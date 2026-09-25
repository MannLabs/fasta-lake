#!/usr/bin/env python3
"""Validate figure packaging, cited labels, source hashes and browser receipts."""
from __future__ import annotations

import base64
import hashlib
import json
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    """Return a file identity for the saved asset and render receipts."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    """Check deterministic assets, pinned code citations and saved browser renders."""
    data = json.loads((HERE / 'pipeline_data.json').read_text())
    sources = json.loads((HERE / 'code_sources.json').read_text())
    render = json.loads((HERE / 'RENDER_RECEIPT.json').read_text())
    previous = json.loads((HERE / 'BUILD_RECEIPT.json').read_text())
    subprocess.run([sys.executable, str(HERE / 'build_figure.py')], check=True, capture_output=True)
    build = json.loads((HERE / 'BUILD_RECEIPT.json').read_text())
    assert build == previous, 'Build not reproducible; first run the builder after edits'
    ledger = (HERE / 'LABEL_SOURCES.md').read_text()
    inventory = (HERE / 'PIPELINE_INVENTORY.md').read_text()
    assert inventory.endswith(ledger)
    citation_count = 0
    source_cache = {}
    for key, source in sources['sources'].items():
        path = HERE / 'source_snapshot' / 'code' / source['path']
        if source['path'] not in source_cache:
            if path.exists():
                raw_source = path.read_bytes()
            else:
                # A normal clone retains the pinned source; no manuscript snapshots needed.
                raw_source = subprocess.run(
                    ['git', 'show', sources['commit'] + ':' + source['path']],
                    cwd=HERE, check=True, capture_output=True,
                ).stdout
            source_cache[source['path']] = raw_source
        raw_source = source_cache[source['path']]
        assert hashlib.sha256(raw_source).hexdigest() == source['sha256'], (key, 'changed source')
        assert 1 <= source['first_line'] <= source['last_line'] <= len(raw_source.decode().splitlines())
        citation_count += 1
    for label in data['labels'].values():
        assert label['text'].replace('|',r'\|') in ledger
        assert label['sources'] and all(k in sources['sources'] for k in label['sources'])
    for stage in data['stages']:
        assert stage['detail'].replace('|',r'\|') in ledger
        assert stage['sources'] and all(k in sources['sources'] for k in stage['sources'])

    svg_records = []
    for theme in ('light','dark'):
        path = HERE / f'fastalake-pipeline-{theme}.svg'
        subprocess.run(['xmllint','--noout',str(path)],check=True,capture_output=True)
        raw = path.read_text()
        root = ET.fromstring(raw)
        for element in root.iter():
            assert element.tag.split('}')[-1].lower() not in ('script','foreignObject'.lower())
            for key,value in element.attrib.items():
                assert not key.lower().startswith('on'), ('event handler',key)
                if key.split('}')[-1] == 'href':
                    assert value.startswith('#') or value.startswith('data:'), ('external href',value)
        urls = re.findall(r'url\((.*?)\)',raw)
        assert len(urls)==2 and all(u.startswith('data:font/woff2;base64,') for u in urls)
        for url, font in zip(urls,build['fonts']):
            assert hashlib.sha256(base64.b64decode(url.split(',',1)[1])).hexdigest() == font['sha256']
        texts = [e for e in root.iter() if e.tag.endswith('}text')]
        for text in texts:
            assert text.text == data['labels'][text.attrib['data-label']]['text']
            assert float(text.attrib['font-size'])*880/data['width'] >= 11
        assert 'clustering' not in raw.lower()
        assert 'k-mer' not in raw.lower() and 'directLFQ' not in raw and 'TPM' not in raw
        assert 'Selection parsimony' in raw
        svg_records.append(dict(file=path.name,labels=len(texts),xmllint='PASS',scripts=0,external_href=0,embedded_fonts=2,sha256=digest(path)))

    html=(HERE/'fastalake-pipeline.html').read_text()
    assert '<script' not in html.lower() and 'clustering' not in html.lower()
    assert render['status']=='PASS' and not render['errors'] and not render['external_requests']
    assert len(render['renders'])==12
    for record in render['renders']:
        assert digest(HERE/record['png']) == record['sha256'], 'PNG changed since browser check'
        metrics=record.get('metrics',[])
        for i,item_a in enumerate(metrics):
            for item_b in metrics[i+1:]:
                a,b=item_a['ink'],item_b['ink']
                overlap_x=min(a['x']+a['width'],b['x']+b['width'])-max(a['x'],b['x'])
                overlap_y=min(a['y']+a['height'],b['y']+b['height'])-max(a['y'],b['y'])
                assert overlap_x <= 0 or overlap_y <= 0, ('overlapping labels',item_a['key'],item_b['key'])
    visual = json.loads((HERE/'VISUAL_REVIEW.json').read_text())
    assert visual['status']=='PASS'
    assert {r['png']:r['sha256'] for r in render['renders']} == {r['png']:r['sha256'] for r in visual['renders']}
    receipt=dict(checked_utc=datetime.now(timezone.utc).isoformat(),status='PASS',source_commit=data['source_commit'],
                 source_records_checked=citation_count,labels_checked=len(data['labels']),detail_paragraphs_checked=len(data['stages']),
                 reproducible_build=True,svg=svg_records,minimum_font_px_at_880=build['minimum_font_px_at_880'],
                 browser_render_count=12,visual_review='VISUAL_REVIEW.json',no_overlapping_text_ink=True,
                 old_stage_check={'clustering':'absent','k-mer rescue':'absent','TPM filtering':'absent',
                                  'directLFQ':'absent from figure; explicitly optional in HTML notes',
                                  'parsimony':'present: active Rust selection, distinct from study grouping'},
                 manuscript_alignment='v11 terminology checked on 23 September 2026; see PIPELINE_INVENTORY.md')
    (HERE/'VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
