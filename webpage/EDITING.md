# BS-Cloth Project Page — Editing Guide

The page is populated from the SIGGRAPH 2026 paper *"Efficient B-Spline Finite
Elements for Cloth Simulation."* Author list, abstract, figure stills, and the
section structure are in place. What remains is the URLs, the YouTube link, and
the videos.

Search `TODO:` in `index.html` to find every editable spot (~20 markers).

---

## 1. Abstract

**File:** `docs/index.html`
**Section:** `<section class="section" id="abstract">`

The three `<p>` tags inside `<div class="content has-text-justified">` are the
three abstract paragraphs. Edit them directly — no special markers, just the
prose itself.

## 2. Resource buttons (arXiv / Paper / Supplemental / Code / Video)

**Search:** `TODO: links`

Five `<span class="link-block">` blocks. Replace each `href="#"` with a real
URL once available. The "Video" button anchors to `#video` (the supplementary
video at the bottom) — already wired, no change needed.

## 3. Canonical / OG URLs

**Search:** `TODO-username`

Replace `TODO-username` everywhere in `<head>` with your actual GitHub
username (or custom domain).

## 4. Supplementary video (YouTube link)

**Search:** `TODO_YOUTUBE_URL`

Single occurrence in the `#video` section near the bottom. Replace with the
full YouTube watch URL (e.g. `https://www.youtube.com/watch?v=XXXXXXXXXXX`).
The teaser image acts as a thumbnail with a red play overlay; clicking opens
YouTube in a new tab.

## 5. Square-cloth comparison videos (12 files)

The "Side-by-Side Video Comparison" table expects videos at
`static/videos/comparison/<scene>_<method>.mp4`. Until they exist, the
per-test poster image (the same still used in the Comparison section above)
is shown in each cell.

| Scene slug | Methods (4 each) |
|---|---|
| `upright_hanging` | `bspline`, `linear`, `sfem`, `bhem` |
| `drape` | `bspline`, `linear`, `sfem`, `bhem` |
| `shear` | `bspline`, `linear`, `sfem`, `bhem` |

Full paths (12 files):
```
static/videos/comparison/upright_hanging_bspline.mp4
static/videos/comparison/upright_hanging_linear.mp4
static/videos/comparison/upright_hanging_sfem.mp4
static/videos/comparison/upright_hanging_bhem.mp4
static/videos/comparison/drape_bspline.mp4
static/videos/comparison/drape_linear.mp4
static/videos/comparison/drape_sfem.mp4
static/videos/comparison/drape_bhem.mp4
static/videos/comparison/shear_bspline.mp4
static/videos/comparison/shear_linear.mp4
static/videos/comparison/shear_sfem.mp4
static/videos/comparison/shear_bhem.mp4
```

To rename slugs/methods, edit the `comparisonScenes` and `comparisonMethods`
arrays in the `<script>` block at the bottom of `index.html`.

## 6. Stress test videos (4 files, carousel)

The Stress Tests section is a horizontal carousel (left/right arrow buttons,
scroll-snap). One card per scene; videos at:

| Path | Scene |
|---|---|
| `static/videos/rotating_sphere.mp4` | Cloth on rotating sphere |
| `static/videos/helicopter.mp4` | Helicopter |
| `static/videos/stripes.mp4` | Stripes sweeping over armadillo |
| `static/videos/cloth_basket.mp4` | Cloth basket (matches teaser) |

## 7. BibTeX

**Search:** `TODO` inside the `<pre><code>@article{` block.

Fill in `volume`, `number`, and `doi` once the publisher assigns them.

## 8. Footer attribution

**Search:** `&copy; 2026 The BS-Cloth authors.`

Replace with your preferred attribution (lab, group name, etc.). Keep the
Nerfies and PAT3D template links — they're a license requirement.

---

## Section structure (current)

| Section ID | Content |
|---|---|
| `#abstract` | Abstract paragraphs from the paper |
| `#comparison-table` | PAT3D-style table — rows = test cases (upright hanging / drape / shear), columns = methods (B-spline / Linear / SFEM / BHEM), cells = videos |
| `#stress` | Horizontal scroll-snap carousel with ‹/› arrow buttons; four scenes (rotating sphere, helicopter, stripes, cloth basket) |
| `#video` | Supplementary video as a clickable YouTube thumbnail near the bottom |
| `#citation` | BibTeX |

All section anchors are wired into the navbar.

## Re-rendering paper figures at higher resolution

The current PNGs were rendered at 120–150 DPI from the paper PDFs:
```bash
cd C:/_dev/BSFEM-Paper/img
pdftoppm -png -r 200 teaser-v2.pdf C:/_dev/BS-Cloth/docs/static/images/teaser
mv C:/_dev/BS-Cloth/docs/static/images/teaser-1.png \
   C:/_dev/BS-Cloth/docs/static/images/teaser.png
```

## Styling

- Project-specific overrides: `docs/stylesheet.css`
- Per-page CSS for cards/tables/carousel/video-link: the `<style>` block in `<head>`
- Bulma framework defaults: `docs/static/css/bulma.min.css` (don't edit Bulma)
