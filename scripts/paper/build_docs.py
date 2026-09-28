"""Build a print-ready PDF source and a valid EPUB 3 from the preprint HTML.

  python build_docs.py <preprint.html> <assets_dir> <out_dir>

The artifact page is the single source of the text. Two derivatives are produced:

  print.html  a standalone document with system fonts, a forced light palette and page rules, ready
              for a headless-browser print. The video becomes a still.
  .epub       EPUB 3.0, XHTML serialised so named entities become literal characters, with the two
              trajectories embedded as animated GIFs. EPUB readers support GIF far more consistently
              than they support video, which is the reason for the format choice.
"""
import html
import os
import re
import shutil
import sys
import uuid
import zipfile
from datetime import datetime, timezone

import lxml.etree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lxml.html
from PIL import Image

import arxiv_css

SRC, ASSETS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
raw = open(SRC, encoding="utf-8").read()

TITLE = re.search(r"<title>(.*?)</title>", raw, re.S).group(1).strip()
body = raw[raw.index('<div class="sheet">'):]
body = body[:body.rindex("</div>") + 6]

AUTHORS = "Joseph Campagna"
UID = "urn:uuid:" + str(uuid.uuid5(uuid.NAMESPACE_URL, "motors-in-the-membrane/preprint/v1"))
NOW = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

BASE_CSS = """
body{font-family:Georgia,"Times New Roman",serif;line-height:1.55;color:#14161a;background:#fff;margin:0}
.sheet{max-width:46em;margin:0 auto;padding:0 1em}
.banner{background:#eceadf;padding:.5em .8em;font-family:"Segoe UI",Helvetica,Arial,sans-serif;
  font-size:.62em;letter-spacing:.1em;text-transform:uppercase;color:#414652;border-bottom:1px solid #d9dad2;margin-bottom:1.4em}
.banner span{margin-right:1.4em;white-space:nowrap}
h1{font-size:1.95em;line-height:1.15;margin:.2em 0 .15em;font-weight:700}
.subtitle{font-size:1.12em;color:#414652;margin:0 0 1em}
.authors{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.92em;margin:0 0 .2em}
.affil,.corr{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.76em;color:#6f7480;margin:0 0 .15em}
.meta{font-family:"Courier New",monospace;font-size:.72em;color:#6f7480;margin:1em 0 0;padding-top:.7em;border-top:1px solid #e9eae3}
.meta span{display:block;margin-bottom:.2em}
.meta b{color:#414652;font-weight:600}
.abstract{background:#f6f6f2;border:1px solid #e9eae3;border-left:3px solid #7a2518;padding:1em 1.2em;margin:1.4em 0 .5em;font-size:.95em}
.abstract h2{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.68em;letter-spacing:.12em;
  text-transform:uppercase;color:#6f7480;margin:0 0 .5em;font-weight:700}
.abstract p{margin:0 0 .6em}
.kw{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.8em;color:#414652;margin:.6em 0 0}
h2.sec{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:1.22em;margin:1.7em 0 .4em;font-weight:600}
h3{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.98em;margin:1.1em 0 .3em;font-weight:600}
p{margin:0 0 .75em}
ul,ol{margin:0 0 .8em;padding-left:1.4em}
li{margin-bottom:.35em}
figure{margin:1.4em 0}
figure img{width:100%;border:1px solid #e9eae3}
figcaption{font-size:.82em;color:#414652;margin-top:.5em;line-height:1.45}
.tablewrap{margin:1.1em 0;overflow-x:auto}
table{border-collapse:collapse;width:100%;font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.78em}
caption{caption-side:top;text-align:left;font-size:.84em;color:#414652;margin-bottom:.5em;line-height:1.45}
th,td{text-align:left;padding:.4em .8em .4em 0;border-bottom:1px solid #e9eae3;vertical-align:top}
thead th{border-bottom:1px solid #d9dad2;font-size:.85em;letter-spacing:.04em;text-transform:uppercase;color:#414652}
td.n,th.n{text-align:right;white-space:nowrap;padding-right:1.2em}
.callout{border:1px solid #d9dad2;border-left:3px solid #1d4e89;background:#f6f6f2;padding:.8em 1em;margin:1.1em 0;font-size:.95em}
.callout p{margin:0}
code{font-family:"Courier New",monospace;font-size:.9em;background:#f3f3ef;padding:0 .2em}
.refs{font-size:.82em;padding-left:1.8em}
.refs li{margin-bottom:.4em;color:#414652}
.foot{font-family:"Segoe UI",Helvetica,Arial,sans-serif;font-size:.76em;color:#6f7480;border-top:1px solid #d9dad2;margin-top:1.8em;padding-top:.8em}
.note{font-size:.82em;color:#6f7480;font-style:italic}
sup.cite{font-size:.7em;color:#1d4e89}
a{color:#1d4e89}
"""

PRINT_CSS = BASE_CSS + arxiv_css.CSS

EPUB_CSS = BASE_CSS + """
html{font-size:100%}
.sheet{padding:0}
figure img{page-break-inside:avoid}
"""


def shrink(src, dst, width):
    im = Image.open(src)
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    im.convert("RGB").save(dst, quality=88, optimize=True)
    return dst


# ---------------------------------------------------------------- assets
img_dir = os.path.join(OUT, "_img")
os.makedirs(img_dir, exist_ok=True)
FIGS = {"fig1-cross-section.png": ("bilayer-cross-section.png", 1500),
        "fig2-tilt-depth.png": ("tilt-and-depth.png", 1500),
        "fig3-torsions.png": ("torsion-validation.png", 1500)}
for out_name, (src_name, w) in FIGS.items():
    src = os.path.join(ASSETS, src_name)
    shrink(src, os.path.join(img_dir, out_name.replace(".png", ".jpg")), w)
GIFS = ["video-s1-closeup.gif", "video-s2-cross-section.gif"]
for g in GIFS:
    shutil.copy(os.path.join(ASSETS, g), os.path.join(img_dir, g))

VIDEO_BLOCK = re.search(r"<figure>\s*<video.*?</figure>", body, re.S).group(0)


HEAD_RE = re.compile(r'<div class="banner">.*?</p>\s*(?=<div class="meta">)', re.S)


def rewrite(markup, for_epub):
    m = markup
    if not for_epub:
        m = HEAD_RE.sub(arxiv_css.TITLE_BLOCK.format(stamp=arxiv_css.STAMP), m, count=1)
    for name in FIGS:
        m = m.replace('src="media/%s"' % name, 'src="%s%s"' % ("images/" if for_epub else "_img/",
                                                               name.replace(".png", ".jpg")))
    if for_epub:
        vid = ('<figure>\n<img src="images/video-s1-closeup.gif" alt="Animated trajectory of the motor"/>\n'
               '<figcaption><b>Supplementary Video S1.</b> MM1 replica 1, camera tracking the ligand centre of '
               'mass, rendered as an animated figure. One second corresponds to 5 ns. Displayed coordinates are '
               'smoothed over 1 ns; analysis uses unsmoothed data.</figcaption>\n</figure>\n'
               '<figure>\n<img src="images/video-s2-cross-section.gif" alt="Animated bilayer cross-section"/>\n'
               '<figcaption><b>Supplementary Video S2.</b> The same trajectory as a cross-section through the '
               'full bilayer, drawn in periodic copies. The scale bar is exact under orthographic projection.'
               '</figcaption>\n</figure>')
    else:
        vid = ('<figure>\n<img src="_img/fig1-cross-section.jpg" alt="Still from the trajectory animation"/>\n'
               '<figcaption><b>Supplementary Videos S1 and S2.</b> Animated trajectories are not reproducible in '
               'print. They are provided as animated figures in the EPUB 3 edition and as MP4 files in the code '
               'repository.</figcaption>\n</figure>')
    return m.replace(VIDEO_BLOCK, vid)


# ---------------------------------------------------------------- PDF source
print_html = ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"/>"
              "<title>%s</title><style>%s</style></head><body>\n%s\n</body></html>"
              % (html.escape(TITLE), PRINT_CSS, rewrite(body, False)))
print_path = os.path.join(OUT, "print.html")
open(print_path, "w", encoding="utf-8").write(print_html)
print("wrote", print_path)

# ---------------------------------------------------------------- EPUB 3
frag = lxml.html.fragment_fromstring(rewrite(body, True))
content = ET.tostring(frag, encoding="unicode", method="xml")
content = re.sub(r'\s(controls|playsinline|muted|loop)(?=[\s>])', '', content)

xhtml = ('<?xml version="1.0" encoding="utf-8"?>\n'
         '<html xmlns="http://www.w3.org/1999/xhtml" xmlns:epub="http://www.idpf.org/2007/ops" lang="en" xml:lang="en">\n'
         '<head><meta charset="utf-8"/><title>%s</title>'
         '<link rel="stylesheet" type="text/css" href="../css/style.css"/></head>\n'
         '<body>\n%s\n</body></html>' % (html.escape(TITLE), content))

nav = ('<?xml version="1.0" encoding="utf-8"?>\n'
       '<html xmlns="http://www.w3.org/1999/xhtml" xmlns:epub="http://www.idpf.org/2007/ops" lang="en" xml:lang="en">\n'
       '<head><meta charset="utf-8"/><title>Contents</title></head><body>\n'
       '<nav epub:type="toc" id="toc"><h1>Contents</h1><ol>\n'
       + "\n".join('<li><a href="text/preprint.xhtml">%s</a></li>' % s for s in
                   ["Abstract", "1. Introduction", "2. Methods", "3. Results", "4. Discussion",
                    "5. Limitations", "6. Conclusions", "Tables", "Supplementary material", "References"])
       + '\n</ol></nav>\n<nav epub:type="landmarks" hidden="hidden"><ol>'
         '<li><a epub:type="bodymatter" href="text/preprint.xhtml">Start</a></li></ol></nav>\n'
         '</body></html>')

items, spine_imgs = [], []
for name in sorted(os.listdir(img_dir)):
    mt = "image/gif" if name.endswith(".gif") else "image/jpeg"
    items.append('<item id="%s" href="images/%s" media-type="%s"/>'
                 % (re.sub(r"\W", "_", name), name, mt))
    spine_imgs.append(name)

opf = ('<?xml version="1.0" encoding="utf-8"?>\n'
       '<package xmlns="http://www.idpf.org/2007/opf" version="3.0" unique-identifier="pub-id" xml:lang="en">\n'
       '<metadata xmlns:dc="http://purl.org/dc/elements/1.1/">\n'
       '<dc:identifier id="pub-id">%s</dc:identifier>\n'
       '<dc:title>%s</dc:title>\n<dc:language>en</dc:language>\n'
       + "".join('<dc:creator>%s</dc:creator>\n' % html.escape(a.strip()) for a in AUTHORS.split(";")) +
       '<dc:date>%s</dc:date>\n<dc:publisher>Preprint, not peer reviewed</dc:publisher>\n'
       '<dc:description>Ten-replica atomistic simulation of light-driven molecular motors in a '
       'bacterial-mimetic bilayer, with a measured sampling floor and a quantum-chemical force-field audit.</dc:description>\n'
       '<meta property="dcterms:modified">%s</meta>\n</metadata>\n'
       '<manifest>\n<item id="nav" href="nav.xhtml" media-type="application/xhtml+xml" properties="nav"/>\n'
       '<item id="main" href="text/preprint.xhtml" media-type="application/xhtml+xml"/>\n'
       '<item id="css" href="css/style.css" media-type="text/css"/>\n%s\n</manifest>\n'
       '<spine>\n<itemref idref="main"/>\n</spine>\n</package>') % (
    UID, html.escape(TITLE), NOW[:10], NOW, "\n".join(items))

epub_path = os.path.join(OUT, "motors-in-the-membrane.epub")
with zipfile.ZipFile(epub_path, "w") as z:
    z.writestr(zipfile.ZipInfo("mimetype"), "application/epub+zip", compress_type=zipfile.ZIP_STORED)
    z.writestr("META-INF/container.xml",
               '<?xml version="1.0" encoding="utf-8"?>\n<container version="1.0" '
               'xmlns="urn:oasis:names:tc:opendocument:xmlns:container"><rootfiles>'
               '<rootfile full-path="EPUB/package.opf" media-type="application/oebps-package+xml"/>'
               '</rootfiles></container>', zipfile.ZIP_DEFLATED)
    z.writestr("EPUB/package.opf", opf, zipfile.ZIP_DEFLATED)
    z.writestr("EPUB/nav.xhtml", nav, zipfile.ZIP_DEFLATED)
    z.writestr("EPUB/text/preprint.xhtml", xhtml, zipfile.ZIP_DEFLATED)
    z.writestr("EPUB/css/style.css", EPUB_CSS, zipfile.ZIP_DEFLATED)
    for name in spine_imgs:
        z.write(os.path.join(img_dir, name), "EPUB/images/" + name,
                zipfile.ZIP_STORED if name.endswith(".gif") else zipfile.ZIP_DEFLATED)
print("wrote %s  (%.1f MB)" % (epub_path, os.path.getsize(epub_path) / 1048576))

# quick structural self-check
with zipfile.ZipFile(epub_path) as z:
    assert z.namelist()[0] == "mimetype", "mimetype must be the first entry"
    assert z.getinfo("mimetype").compress_type == zipfile.ZIP_STORED, "mimetype must be stored"
    ET.fromstring(z.read("EPUB/text/preprint.xhtml"))
    ET.fromstring(z.read("EPUB/package.opf"))
    ET.fromstring(z.read("EPUB/nav.xhtml"))
    print("  mimetype stored and first; XHTML, OPF and nav are well-formed XML")
    print("  animated figures:", ", ".join(n for n in z.namelist() if n.endswith(".gif")))
