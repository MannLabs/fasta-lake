/* Progressive enhancement: the first frame and all prose remain without JS. */
(() => {
  function initialise() {
    document.querySelectorAll('.fl-explainer:not([data-ready])').forEach(async box => {
      box.dataset.ready = 'true';
      const picture = box.querySelector('img');
      if (!picture) return;
      try {
        const base = new URL('.', picture.src.startsWith('data:') ? document.baseURI : picture.src);
        const embedded = box.querySelector('script[type="application/json"]');
        let frames;
        if (embedded) frames = JSON.parse(embedded.textContent);
        else {
          const response = await fetch(new URL(box.dataset.story, base));
          if (!response.ok) throw new Error('Storyboard unavailable');
          frames = await response.json();
        }
        if (!frames.length) return;
        let index = 0, timer = null, objectURL = null;
        const controls = document.createElement('div');
        controls.className = 'fl-controls';
        const caption = document.createElement('p');
        caption.className = 'fl-caption';
        caption.setAttribute('aria-live', 'polite');
        const progress = document.createElement('span');
        const fullFrame = document.createElement('a');
        fullFrame.textContent = 'Open full-size frame';
        fullFrame.target = '_blank';
        fullFrame.rel = 'noopener';
        const button = (label, action) => {
          const node = document.createElement('button');
          node.type = 'button'; node.textContent = label;
          node.addEventListener('click', action); controls.append(node);
          return node;
        };
        const stop = () => {
          clearInterval(timer); timer = null; play.textContent = 'Play';
          play.setAttribute('aria-pressed', 'false');
        };
        const show = () => {
          picture.src = new URL(frames[index].image, base).href;
          if (objectURL) URL.revokeObjectURL(objectURL);
          objectURL = null;
          if (picture.src.startsWith('data:image/svg+xml;base64,')) {
            const bytes = Uint8Array.from(atob(picture.src.split(',')[1]), char => char.charCodeAt(0));
            objectURL = URL.createObjectURL(new Blob([bytes], {type: 'image/svg+xml'}));
          }
          fullFrame.href = objectURL || picture.src;
          picture.alt = frames[index].caption;
          caption.textContent = frames[index].caption;
          progress.textContent = `Step ${index + 1} of ${frames.length}`;
          previous.disabled = index === 0;
          next.disabled = index === frames.length - 1;
        };
        const previous = button('Previous', () => { stop(); index--; show(); });
        const play = button('Play', () => {
          if (timer) { stop(); return; }
          if (index === frames.length - 1) index = 0;
          show(); play.textContent = 'Pause'; play.setAttribute('aria-pressed', 'true');
          timer = setInterval(() => {
            if (!box.isConnected || document.hidden) { stop(); return; }
            index++; show();
            if (index === frames.length - 1) stop();
          }, 4500);
        });
        play.setAttribute('aria-pressed', 'false');
        const next = button('Next', () => { stop(); index++; show(); });
        button('Restart', () => { stop(); index = 0; show(); });
        controls.append(progress);
        box.append(controls, caption, fullFrame);
        // No autoplay, including for visitors preferring reduced motion.
        show();
      } catch (_) {
        // Keep the static image and written steps usable if offline/fetch fails.
      }
    });
  }
  if (typeof document$ !== 'undefined') document$.subscribe(initialise);
  else if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initialise);
  else initialise();
})();
