#!/usr/bin/env python3
"""Draw SNFG diagrams from canonical topology strings (topoSeq), in the style of output/snfg_by_g (horizontal,
reducing end at the right, no linkage labels). Used for the species that have no GlyTouCan identifier.
Usage: draw_toposeq_snfg.py manual_structures.csv out_dir   (CSV columns: file,toposeq)"""
import csv, sys, pathlib
COL = {"GlcNAc": "#0385AE", "GalNAc": "#FCC326", "Man": "#058F60", "Gal": "#FCC326", "Glc": "#0385AE", "Neu5Ac": "#A15989",
       "Neu5Gc": "#91D3E3", "Fuc": "#C23537", "Xyl": "#F47F2A", "Hex": "#FFFFFF", "HexNAc": "#FFFFFF"}
SHAPE = {"GlcNAc": "square", "GalNAc": "square", "Man": "circle", "Gal": "circle", "Glc": "circle", "Neu5Ac": "diamond",
         "Neu5Gc": "diamond", "Fuc": "triangle", "Xyl": "star", "Hex": "circle", "HexNAc": "square"}
U = 100  # spacing between residues

def parse(s):
    """'(A(B)(C(D)))' -> ('A', [('B', []), ('C', [('D', [])])])"""
    pos = 0
    def node():
        nonlocal pos
        assert s[pos] == "("; pos += 1
        j = pos
        while s[j] not in "()": j += 1
        name = s[pos:j]; pos = j; kids = []
        while s[pos] == "(":
            kids.append(node())
        assert s[pos] == ")"; pos += 1
        return (name, kids)
    return node()

def layout(tree):
    """x = -depth*U (root at 0, children to the left); y from leaf order, parents centred on their children."""
    nodes, edges = [], []
    leaf = [0]
    def place(n, depth, parent):
        name, kids = n
        if name == "Fuc" and parent is not None and depth == 1:        # core fucose hangs below the reducing-end GlcNAc
            idx = len(nodes); nodes.append([name, parent[1], parent[2] + U]); edges.append((parent[0], idx)); return idx
        if not kids:
            y = leaf[0] * U; leaf[0] += 1
            idx = len(nodes); nodes.append([name, -depth * U, y])
            if parent is not None: edges.append((parent[0], idx))
            return idx
        idx = len(nodes); nodes.append([name, -depth * U, 0.0]); me = (idx, -depth * U, 0.0)
        if parent is not None: edges.append((parent[0], idx))
        ys = []
        for k in kids:
            if k[0] == "Fuc" and depth == 0:
                ys.append(None); place(k, depth + 1, me); continue
            cid = place(k, depth + 1, me); ys.append(nodes[cid][2])
        ys = [y for y in ys if y is not None]
        nodes[idx][2] = sum(ys) / len(ys) if ys else 0.0
        return idx
    place(tree, 0, None)
    # Fuc placed below the root: fix its y after the root's y is known
    root_y = nodes[0][2]
    for a, b in edges:
        if a == 0 and nodes[b][0] == "Fuc" and nodes[b][1] == 0: nodes[b][2] = root_y + U
    return nodes, edges

def svg(nodes, edges):
    xs = [n[1] for n in nodes]; ys = [n[2] for n in nodes]
    x0, x1, y0, y1 = min(xs) - 50, max(xs) + 50, min(ys) - 50, max(ys) + 50
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{x1-x0:.0f}" height="{y1-y0:.0f}" viewBox="{x0:.0f} {y0:.0f} {x1-x0:.0f} {y1-y0:.0f}">']
    for a, b in edges:
        out.append(f'<path d="M{nodes[a][1]},{nodes[a][2]} L{nodes[b][1]},{nodes[b][2]}" stroke-width="4" stroke="#000"/>')
    for name, x, y in nodes:
        c, sh = COL.get(name, "#FFFFFF"), SHAPE.get(name, "circle"); st = 'stroke-width="2" stroke="#000"'
        if sh == "square": out.append(f'<rect x="{x-25}" y="{y-25}" width="50" height="50" fill="{c}" {st}/>')
        elif sh == "circle": out.append(f'<circle cx="{x}" cy="{y}" r="25" fill="{c}" {st}/>')
        elif sh == "diamond": out.append(f'<polygon points="{x},{y-28} {x+28},{y} {x},{y+28} {x-28},{y}" fill="{c}" {st}/>')
        elif sh == "triangle": out.append(f'<polygon points="{x},{y-27} {x+27},{y+22} {x-27},{y+22}" fill="{c}" {st}/>')
        else: out.append(f'<circle cx="{x}" cy="{y}" r="25" fill="{c}" {st}/>')
    out.append("</svg>"); return "\n".join(out)

if __name__ == "__main__":
    src, out = sys.argv[1], pathlib.Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    n = 0
    for row in csv.DictReader(open(src, encoding="utf-8")):
        nodes, edges = layout(parse(row["toposeq"]))
        (out / row["file"]).write_text(svg(nodes, edges)); n += 1
    print("drew", n, "structures into", out)
