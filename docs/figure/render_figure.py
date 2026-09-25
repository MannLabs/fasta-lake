#!/usr/bin/env python3
"""Render the SVG/HTML matrix and audit browser geometry without network access."""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent


def main():
    """Render both themes at two sizes and save geometry and asset identities."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromium", help="Optional existing Chromium executable")
    args = parser.parse_args()
    out = HERE / "renders"
    out.mkdir(exist_ok=True)
    bounds = json.loads((HERE / "LAYOUT_BOUNDS.json").read_text())
    records, requests, errors = [], [], []
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=args.chromium,
                                    args=["--no-sandbox", "--disable-dev-shm-usage"])
        for width in (880, 1600):
            for scheme in ("light", "dark"):
                context = browser.new_context(viewport={"width": width, "height": math.ceil(width * .9)},
                                              device_scale_factor=1, color_scheme=scheme)
                def block(route):
                    """Record and block external requests during the offline render check."""
                    requests.append(route.request.url)
                    route.abort()
                context.route("http://**/*", block)
                context.route("https://**/*", block)
                page = context.new_page()
                page.on("pageerror", lambda error: errors.append(str(error)))
                for theme in ("light", "dark"):
                    filename = f"fastalake-pipeline-{theme}.svg"
                    svg = (HERE / filename).read_text()
                    # GitHub embeds SVG as an image, where font loading has stricter rules.
                    uri = "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()
                    page.set_content('<html><body style="margin:0"><img style="display:block;width:100%;height:auto" '
                                     f'src="{uri}"></body></html>')
                    page.locator("img").evaluate("img => img.decode()")
                    path = out / f"svg-{theme}-{width}-browser-{scheme}.png"
                    page.screenshot(path=str(path), full_page=True)
                    # Inline rendering allows font and geometry checks. It must match img mode.
                    page.set_content('<html><body style="margin:0">' + svg + '</body></html>')
                    page.evaluate("document.querySelector('svg').style.cssText='display:block;width:100%;height:auto'")
                    page.evaluate("document.fonts.ready")
                    font_status = page.evaluate("""() => [...document.fonts].map(f => ({family:f.family,status:f.status}))""")
                    metrics = page.evaluate("""() => [...document.querySelectorAll('text[data-label]')].map(t => {
                        const b=t.getBBox(),ctx=document.createElement('canvas').getContext('2d');
                        ctx.font=`${t.getAttribute('font-weight')} ${t.getAttribute('font-size')}px "${t.getAttribute('font-family')}"`;
                        const m=ctx.measureText(t.textContent),x=Number(t.getAttribute('x')),y=Number(t.getAttribute('y'));
                        return {key:t.dataset.label, x:b.x,y:b.y,width:b.width,height:b.height,
                          ink:{x:x-m.actualBoundingBoxLeft,y:y-m.actualBoundingBoxAscent,
                            width:m.actualBoundingBoxLeft+m.actualBoundingBoxRight,
                            height:m.actualBoundingBoxAscent+m.actualBoundingBoxDescent}};
                    })""")
                    inline = page.screenshot(full_page=True)
                    # Raster equality confirms font subsets work in both embedding contexts.
                    same_pixels = inline == path.read_bytes()
                    record = dict(asset=filename,width=width,browser_scheme=scheme,png=str(path.relative_to(HERE)),
                                  sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                  image_and_inline_identical=same_pixels,fonts=font_status,metrics=metrics)
                    for item in metrics:
                        spec = next(b for b in bounds if b['theme']==theme and b['key']==item['key'])
                        if item['width'] > spec['max_width'] + .5:
                            errors.append(f"{filename}/{width}/{item['key']}: text {item['width']:.2f} exceeds {spec['max_width']}")
                        if item['x'] < 0 or item['y'] < 0 or item['x']+item['width'] > 1200 or item['y']+item['height'] > 1080:
                            errors.append(f"{filename}/{item['key']}: outside canvas")
                    if any(f['status'] != 'loaded' for f in font_status) or len(font_status) != 2:
                        errors.append(f"{filename}: embedded fonts failed to load")
                    if not same_pixels:
                        errors.append(f"{filename}/{width}: img and inline render differ")
                    records.append(record)
                page.goto((HERE / "fastalake-pipeline.html").as_uri())
                page.evaluate("document.fonts.ready")
                visible = page.locator(".theme:visible")
                if visible.count() != 1 or scheme not in visible.get_attribute('class'):
                    errors.append(f"HTML {scheme}: theme switching failed")
                path = out / f"html-{scheme}-{width}.png"
                page.screenshot(path=str(path), full_page=True)
                records.append(dict(asset="fastalake-pipeline.html",width=width,browser_scheme=scheme,
                                    png=str(path.relative_to(HERE)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                    page_width=page.evaluate("document.documentElement.scrollWidth")))
                if page.evaluate("document.documentElement.scrollWidth") != width:
                    errors.append(f"HTML {width}: horizontal overflow")
                if width == 1600 and scheme == "light":
                    # A vector PDF for easy author review; CSS hides below-figure notes.
                    page.pdf(path=str(HERE / "fastalake-pipeline.pdf"), width="1200px", height="1080px",
                             print_background=True, margin={k:"0px" for k in ('top','right','bottom','left')})
                context.close()
        version = browser.version
        browser.close()
    if requests:
        errors.append("External requests attempted")
    receipt = dict(browser=version,renders=records,external_requests=requests,errors=errors,
                   status="PASS" if not errors else "FAIL",visual_inspection="pending; separate human-visible review required")
    (HERE / "RENDER_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k != 'renders'},indent=2))
    print(f"Rendered {len(records)} PNGs")
    raise SystemExit(bool(errors))


if __name__ == "__main__":
    main()
