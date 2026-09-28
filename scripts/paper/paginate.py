"""Stamp centred page numbers onto the rendered PDF, the way an article-class document carries them.

  python paginate.py <in.pdf> <chrome.exe> <work_dir>

Chrome's print path cannot fill CSS paged-media margin boxes, and its own header/footer template is
not reachable from the command line. So the numbers are rendered as a second, otherwise empty PDF of
the same page count and merged underneath. No external typesetting dependency.
"""
import os
import subprocess
import sys

from pypdf import PdfReader, PdfWriter

src, chrome, work = sys.argv[1], sys.argv[2], sys.argv[3]
n = len(PdfReader(src).pages)

pages = "\n".join(
    '<div class="pg"><span>%d</span></div>' % i for i in range(1, n + 1))
html = """<!doctype html><html><head><meta charset="utf-8"/><style>
@page { size: letter; margin: 0; }
html,body { margin:0; padding:0; }
.pg { width:8.5in; height:11in; position:relative; break-after:page; page-break-after:always; }
.pg:last-child { break-after:auto; page-break-after:auto; }
.pg span { position:absolute; left:0; right:0; bottom:0.62in; text-align:center;
           font-family:"Nimbus Roman","Times New Roman",Times,serif; font-size:10pt; color:#000; }
</style></head><body>%s</body></html>""" % pages

overlay_html = os.path.join(work, "_pagenums.html")
overlay_pdf = os.path.join(work, "_pagenums.pdf")
open(overlay_html, "w", encoding="utf-8").write(html)
subprocess.run([chrome, "--headless=new", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
                "--virtual-time-budget=15000", "--print-to-pdf=" + overlay_pdf,
                "file:///" + overlay_html.replace("\\", "/")],
               capture_output=True, check=False)

nums = PdfReader(overlay_pdf)
assert len(nums.pages) == n, "overlay has %d pages, document has %d" % (len(nums.pages), n)
w = PdfWriter()
for page, num in zip(PdfReader(src).pages, nums.pages):
    page.merge_page(num)
    w.add_page(page)
with open(src, "wb") as f:
    w.write(f)
os.remove(overlay_html)
os.remove(overlay_pdf)
print("stamped %d page numbers into %s" % (n, os.path.basename(src)))
