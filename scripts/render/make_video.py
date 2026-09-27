"""Encode a rendered frame sequence to MP4 with a burned-in legend, clock and scale bar.

  python make_video.py <frames_dir> <side|closeup> <out.mp4> [--crf 19] [--width 0] [--ortho 100.8]

The swatch colours are the colours as they appear after Blender's AgX view transform, not the raw
material values, so the key matches what the viewer actually sees.
"""
import os
import subprocess
import sys

FFMPEG = ("C:/Users/ltjjp/AppData/Local/Microsoft/WinGet/Packages/"
          "Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe/ffmpeg-8.0.1-full_build/bin/ffmpeg.exe")
FONT = "C\\:/Windows/Fonts/consola.ttf"
if not os.path.exists("C:/Windows/Fonts/consola.ttf"):
    FONT = "C\\:/Windows/Fonts/arial.ttf"

frames, view, out = sys.argv[1:4]
args = sys.argv[4:]


def opt(name, default):
    return args[args.index(name) + 1] if name in args else default


crf = opt("--crf", "19")
width = int(opt("--width", "0"))          # 0 keeps 1920
ortho = float(opt("--ortho", "100.8"))    # side view width in Angstrom, for the scale bar
W, H = 1920, 1080
PS_NS = 0.2                               # 200 ps per frame

LEGEND = [
    ("0xE8A860", "motor (MM1)"),
    ("0xF0705E", "C=C axle"),
    ("0x9FE3AC", "ammonium N+"),
    ("0xA9DEEA", "phosphate P"),
    ("0xA9B5C6", "lipid tails"),
]
SUB = "POPE:POPG 1:3 bilayer  .  Amber Lipid21 + GAFF2  .  coordinates smoothed over 1 ns"


def esc(t):
    """drawtext needs colons escaped even inside a quoted value."""
    return t.replace("\\", "\\\\").replace(":", "\:").replace("'", "\'")

f = []
clock = ("drawtext=fontfile='%s':text='MM1 replica 1    t = "
         "%%{eif\\:trunc(n/5)\\:d}.%%{eif\\:mod(n*2,10)\\:d} ns'"
         ":x=48:y=40:fontsize=42:fontcolor=white@0.95"
         ":shadowcolor=black@0.6:shadowx=2:shadowy=2") % FONT
f.append(clock)
f.append(("drawtext=fontfile='%s':text='%s':x=48:y=96:fontsize=25:fontcolor=white@0.62"
          ":shadowcolor=black@0.6:shadowx=1:shadowy=1") % (FONT, esc(SUB)))

x = 48
for colour, label in LEGEND:
    f.append("drawbox=x=%d:y=%d:w=22:h=22:color=%s@1:t=fill" % (x, H - 56, colour))
    f.append(("drawtext=fontfile='%s':text='%s':x=%d:y=%d:fontsize=26:fontcolor=white@0.85"
              ":shadowcolor=black@0.6:shadowx=1:shadowy=1") % (FONT, esc(label), x + 32, H - 58))
    x += 34 + int(15.5 * len(label)) + 26

if view == "side":
    px = int(round(10.0 * W / ortho))      # 1 nm in pixels, exact for an orthographic camera
    bx = W - 60 - px
    f.append("drawbox=x=%d:y=%d:w=%d:h=6:color=white@0.9:t=fill" % (bx, H - 52, px))
    f.append("drawbox=x=%d:y=%d:w=3:h=18:color=white@0.9:t=fill" % (bx, H - 58))
    f.append("drawbox=x=%d:y=%d:w=3:h=18:color=white@0.9:t=fill" % (bx + px - 3, H - 58))
    f.append(("drawtext=fontfile='%s':text='1 nm':x=%d:y=%d:fontsize=26:fontcolor=white@0.9"
              ":shadowcolor=black@0.6:shadowx=1:shadowy=1") % (FONT, bx + px // 2 - 26, H - 92))

if width:
    f.append("scale=%d:-2" % width)
f.append("format=yuv420p")

cmd = [FFMPEG, "-y", "-loglevel", "error", "-framerate", "25", "-start_number", "0",
       "-i", os.path.join(frames, "frame_%04d.png"), "-vf", ",".join(f),
       "-c:v", "libx264", "-crf", crf, "-preset", "slow", "-movflags", "+faststart", out]
print("encoding %s -> %s" % (view, os.path.basename(out)))
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode:
    print(r.stderr[-1500:])
    sys.exit(1)
print("  %.1f MB" % (os.path.getsize(out) / 1048576))
