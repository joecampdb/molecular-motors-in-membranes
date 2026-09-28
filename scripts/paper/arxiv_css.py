"""arXiv/LaTeX-article presentation for the print edition, injected by build_docs.py."""

# Deliberately a placeholder. A well-formed arXiv identifier on a document that has not been
# deposited would read as a real record; the X's keep the layout honest.
STAMP = "arXiv:XXXX.XXXXX [q-bio.BM]  27 Sep 2026"

TITLE_BLOCK = """
<div class="arxiv-stamp">{stamp}</div>
<div class="titleblock">
  <h1 class="ttl">Motors in the Membrane</h1>
  <p class="ttl-sub">Ten-replica atomistic simulation of light-driven molecular motors in a
  bacterial-mimetic bilayer, with a measured sampling floor and a quantum-chemical audit of the
  ligand force field</p>
  <p class="au">Joseph Campagna, PSM</p>
  <p class="aff"><i>Freelance Research</i></p>
  <p class="dt">27 September 2026</p>
</div>
"""

CSS = r"""
/* the side inset is body padding rather than page margin, so the arXiv stamp can sit
   beside the text block and still fall inside Chrome's printable page box */
@page { size: letter; margin: 1in 0.55in; }
html { font-size: 10pt; }
body { font-family: "Nimbus Roman", "Times New Roman", Times, serif; line-height: 1.24;
       color: #000; background: #fff; margin: 0; padding: 0 0.5in; text-align: justify;
       -webkit-font-smoothing: antialiased; position: relative; }
.sheet { max-width: none; padding: 0; margin: 0; }

/* the stamp arXiv prints down the left edge of page one */
.arxiv-stamp { position: absolute; left: 0.04in; top: 7.15in; transform: rotate(-90deg);
               transform-origin: 0 0; white-space: nowrap; font-size: 8.2pt; letter-spacing: .01em; }

/* title block, centred, LaTeX article proportions */
.titleblock { text-align: center; margin: 0 0 1.6em; }
h1.ttl { font-size: 17.5pt; font-weight: 700; line-height: 1.16; margin: 0 0 .35em;
         letter-spacing: -.005em; text-align: center; }
.ttl-sub { font-size: 11pt; line-height: 1.3; margin: 0 auto .95em; max-width: 34em;
           text-align: center; font-style: italic; color: #000; }
.au { font-size: 11.5pt; margin: 0 0 .12em; text-align: center; }
.aff { font-size: 10pt; margin: 0 0 .45em; text-align: center; }
.dt { font-size: 10pt; margin: 0; text-align: center; }

/* abstract: narrower measure, centred heading, as article class sets it */
.abstract { background: none; border: none; padding: 0; margin: 0 2.6em 1.4em; font-size: 9.4pt;
            line-height: 1.26; }
.abstract h2 { font-family: inherit; font-size: 10pt; font-weight: 700; letter-spacing: 0;
               text-transform: none; text-align: center; color: #000; margin: 0 0 .5em; }
.abstract p { margin: 0 0 .5em; text-indent: 0; }
.kw { font-family: inherit; font-size: 9.2pt; color: #000; margin: .7em 0 0; text-indent: 0; }

/* sections */
h2.sec { font-family: inherit; font-size: 12pt; font-weight: 700; margin: 1.35em 0 .35em;
         text-align: left; }
h3 { font-family: inherit; font-size: 10.5pt; font-weight: 700; margin: .95em 0 .25em;
     text-align: left; }
h4 { font-family: inherit; font-size: 10pt; font-weight: 700; font-style: italic; margin: .7em 0 .2em; }

/* LaTeX paragraphing: first line indented, no vertical gap, no indent after a heading */
p { margin: 0; text-indent: 1.6em; }
h1 + p, h2 + p, h3 + p, h4 + p, .callout p, figcaption, caption, li p,
.abstract p, .titleblock p, .meta span, .foot, .note { text-indent: 0; }
ul, ol { margin: .35em 0 .45em; padding-left: 1.7em; }
li { margin-bottom: .18em; text-align: justify; }

/* the web metadata strip becomes a compact preprint header line */
.meta { font-family: inherit; font-size: 8.6pt; color: #000; border: none; border-top: .5pt solid #000;
        border-bottom: .5pt solid #000; padding: .4em 0; margin: 0 0 1.3em; text-align: left;
        display: block; }
.meta span { display: block; margin: 0 0 .1em; }
.meta b { font-weight: 700; }
.banner { display: none; }

/* figures */
figure { margin: 1.1em 0; text-align: center; break-inside: avoid; page-break-inside: avoid; }
figure img { width: 100%; max-width: 100%; max-height: 112mm; object-fit: contain; border: none; }
figcaption { font-size: 8.8pt; line-height: 1.26; margin-top: .45em; text-align: justify;
             color: #000; }

/* booktabs rules: no verticals, no zebra */
.tablewrap { margin: 1.1em 0; break-inside: avoid; page-break-inside: avoid; }
table { border-collapse: collapse; width: 100%; font-family: inherit; font-size: 8.8pt;
        border-top: 1pt solid #000; border-bottom: 1pt solid #000; }
caption { caption-side: top; text-align: justify; font-size: 8.8pt; color: #000; margin-bottom: .45em;
          line-height: 1.26; }
thead th { border-bottom: .5pt solid #000; font-family: inherit; font-size: 8.8pt; font-weight: 700;
           letter-spacing: 0; text-transform: none; color: #000; padding: .3em .7em .3em 0; }
td, th { border: none; padding: .26em .7em .26em 0; vertical-align: top; }
td.n, th.n { text-align: right; padding-right: 1.1em; }
tbody tr:last-child td { border-bottom: none; }

/* callouts become plain indented remarks, as a LaTeX quote environment */
.callout { border: none; border-left: 1.5pt solid #000; background: none; padding: 0 0 0 1em;
           margin: .7em 0 .7em 1.4em; font-size: 9.4pt; }
.callout p { text-indent: 0; }

.refs { font-size: 9pt; padding-left: 2em; }
.refs li { margin-bottom: .22em; color: #000; text-align: justify; }
.note { font-size: 8.8pt; color: #000; font-style: italic; }
.foot { font-family: inherit; font-size: 8.6pt; color: #000; border-top: .5pt solid #000;
        margin-top: 1.4em; padding-top: .5em; }
sup.cite { font-size: 7.2pt; color: #000; vertical-align: super; }
code { font-family: "Nimbus Mono PS", "Courier New", monospace; font-size: .88em; background: none;
       padding: 0; }
a { color: #000; text-decoration: none; }
h2.sec, h3, h4 { break-after: avoid; page-break-after: avoid; }
"""
