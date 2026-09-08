"""Render immutable supplied coordinates; ray-traced static views, one CPU thread."""
import argparse, hashlib, json, pathlib
import numpy as np
import pymol2
from PIL import Image, ImageDraw, ImageFont
P=pathlib.Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--out',default=str(P));args=a.parse_args();out=pathlib.Path(args.out);fig=out/'delivery/figures';fig.mkdir(parents=True,exist_ok=True)
source=P/'inputs/structures'; paths=[source/'1h1q.pdb',source/'vendor_pose_priority.sdf']; before={str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
receipt={'source_hashes':before,'display_policy':'All input atoms retained; cartoon focuses on author chain A, near-site residues and waters shown. Other objects/atoms hidden only for display. No coordinate fitting, repair, preparation or dynamics.','views':{}}
with pymol2.PyMOL() as pm:
 c=pm.cmd;c.set('max_threads',1);c.set('ray_shadows',0);c.set('antialias',2);c.set('ray_trace_mode',1);c.set('orthoscopic',1);c.set('specular',0.15);c.bg_color('white')
 c.load(str(paths[0]),'complex');c.load(str(paths[1]),'priority');c.select('lig','complex and chain A and resn 2A6');c.select('site','byres (complex and chain A and polymer within 4.5 of lig)');c.select('wat','complex and resn HOH within 3.5 of lig')
 xyz0=c.get_coords('all').copy();c.hide('everything','all');c.show('cartoon','complex and chain A');c.color('gray80','complex');c.set('cartoon_transparency',0.6);c.show('sticks','site');c.set('stick_radius',0.15,'site');c.color('gray65','site');c.show('sticks','lig');c.set('stick_radius',0.22,'lig');c.color('forest','lig and elem C');c.color('blue','lig and elem N');c.color('red','lig and elem O');c.show('spheres','wat');c.set('sphere_scale',0.20,'wat');c.color('red','wat');c.color('blue','site and elem N');c.color('red','site and elem O')
 c.set('label_color','black');c.set('label_size',20);c.set('label_outline_color','white');c.label('complex and chain A and resi 81+83 and name CA','resn+resi');c.label('wat','"Water "+resi');
 for name,left,right in [('contact81','lig and name N9','complex and chain A and resi 81 and name O'),('contact83a','lig and name N3','complex and chain A and resi 83 and name N'),('contact83b','lig and name N2','complex and chain A and resi 83 and name O')]:c.distance(name,left,right)
 c.hide('cartoon','all');c.hide('sticks','site');c.show('sticks','complex and chain A and resi 81+83');c.set('dash_color','black');c.set('dash_radius',0.035);c.set('label_size',17,'contact*');c.orient('lig');c.turn('x',25);c.turn('z',-20);c.zoom('lig',7)
 receipt['views']['site']=list(c.get_view());c.png(str(fig/'structure_site.png'),width=1300,height=1000,dpi=150,ray=1)
 c.show('cartoon','complex and chain A');c.hide('sticks','complex and chain A and resi 81+83');c.hide('labels','all');c.disable('contact*');c.hide('sticks','site');c.show('sticks','priority');c.color('magenta','priority and elem C');c.color('blue','priority and elem N');c.color('red','priority and elem O');c.set('stick_radius',0.21,'priority');c.set('cartoon_transparency',0.3);c.orient('lig or priority');c.turn('x',30);c.turn('y',20);c.zoom('lig or priority',7);receipt['views']['frame']=list(c.get_view());c.png(str(fig/'structure_frame.png'),width=1300,height=1000,dpi=150,ray=1)
 assert np.array_equal(xyz0,c.get_coords('all'));receipt['pymol_version']=c.get_version()[0];receipt['coordinates_unchanged']=True
 c.save(str(fig/'structure_scene.pse'))
font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'; ftitle=ImageFont.truetype(font,31);ftext=ImageFont.truetype(font,24)
canvas=Image.new('RGB',(2600,1225),'white');canvas.paste(Image.open(fig/'structure_site.png'),(0,105));canvas.paste(Image.open(fig/'structure_frame.png'),(1300,105));draw=ImageDraw.Draw(canvas)
for x,title,sub in [(25,'1H1Q: observed NU6094 / 2A6 site','Green ligand; gray CDK2; red water. Distances in Å.'),(1325,'Supplied priority export: frame failure','Green = deposited; magenta = unchanged priority export.')]:draw.text((x,20),title,font=ftitle,fill='black');draw.text((x,65),sub,font=ftext,fill='#444444')
draw.text((25,1110),'Polar proximity supports a hinge-binding interpretation; it does not establish energetic importance.',font=ftext,fill='black');draw.text((25,1150),'Priority translated 8.775 Å. An independent ligand fit gives 0.000 Å but conceals the site error.',font=ftext,fill='black');draw.text((25,1190),'Fixed experimental/derived coordinates, no simulation. No pose transfer to other candidates is established.',font=ftext,fill='black');canvas.save(fig/'structure_comparison.png')
assert before=={str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
receipt['outputs']={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in fig.glob('structure*.png')};receipt['font']=font;receipt['font_sha256']=hashlib.sha256(pathlib.Path(font).read_bytes()).hexdigest();(out/'analysis/structure_render_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('Rendered and coordinate/hash checks passed')
