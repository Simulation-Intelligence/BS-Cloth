# BS-Cloth — Project Page

Nerfies-style academic project website, adapted from
[PAT3D](https://github.com/Simulation-Intelligence/PAT3D) /
[Nerfies](https://github.com/nerfies/nerfies.github.io).
No build step — pure static HTML/CSS/JS.

## Local preview

```bash
cd docs
python -m http.server 8000
```
Then open `http://localhost:8000/` in your browser.

On Windows you can also double-click `serve.bat`.

To iterate: edit `index.html` or `stylesheet.css`, then **hard-refresh**
the browser (`Ctrl + Shift + R` / `Cmd + Shift + R`).

## Editing guide

See [EDITING.md](EDITING.md) for a map of every `TODO:` placeholder and where
to drop your media files.

## Publish to GitHub Pages

1. Push this repo to GitHub.
2. Go to repo **Settings → Pages**.
3. Set **Source = Deploy from a branch**, **Branch = main**, **Folder = /docs**.
4. Your page will be at `https://<username>.github.io/BS-Cloth/` in ~60 seconds.

## Troubleshooting

| Symptom | Fix |
|---|---|
| Blank page / nothing loads | Open DevTools → Console + Network. Look for 404s or JS errors. |
| Styles broken on GitHub Pages but fine locally | Asset paths must be **relative** (`static/…`, not `/static/…`). Check for leading slashes. |
| Video doesn't play | Confirm the `.mp4` file exists at the path in `<source src="">`. Check file size isn't 0 bytes. |
| Old version shows after editing | Hard-refresh (`Ctrl+Shift+R`). Or in DevTools → Network, tick "Disable cache" while DevTools is open. |
| FontAwesome icons appear as boxes | The JS bundle (`static/js/fontawesome.all.min.js`) must load. Check Network tab for a 404. |
