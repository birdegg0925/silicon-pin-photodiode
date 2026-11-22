"""Generate SVG figures from the saved scalar summary tables."""
from pathlib import Path
import argparse,csv,html,math
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=ROOT/'pin');ap.add_argument('--out',type=Path,default=ROOT/'figures');a=ap.parse_args()
a.out.mkdir(parents=True,exist_ok=True)
rows=[r for r in csv.DictReader((a.data/'receiver_summary.csv').open()) if r['case_id'].startswith('d35.6_')]
xs=[float(r['case_id'].split('_v')[1]) for r in rows]
def plot(name,title,ylab,series,limit):
    x=lambda v:90+(v-1)/7*570
    y=lambda v:350-v/limit*270
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="450" viewBox="0 0 760 450">','<rect width="760" height="450" fill="white"/>','<g font-family="Arial,sans-serif" fill="#222">',f'<text x="380" y="32" text-anchor="middle" font-size="20">{html.escape(title)}</text>']
    for j in range(6):
        value=limit*j/5; yy=y(value)
        parts += [f'<line x1="90" y1="{yy}" x2="660" y2="{yy}" stroke="#ddd"/>',f'<text x="78" y="{yy+5}" text-anchor="end" font-size="13">{value:g}</text>']
    parts += ['<path d="M90 80V350H660" fill="none" stroke="#333"/>',f'<text x="25" y="215" text-anchor="middle" transform="rotate(-90 25 215)">{html.escape(ylab)}</text>','<text x="375" y="396" text-anchor="middle">Reverse bias (V)</text>']
    for v in xs:parts.append(f'<text x="{x(v)}" y="373" text-anchor="middle">{v:g}</text>')
    for k,(label,values,color) in enumerate(series):
        points=' '.join(f'{x(v)},{y(z)}' for v,z in zip(xs,values))
        parts.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.5"/>')
        for v,z in zip(xs,values):parts.append(f'<circle cx="{x(v)}" cy="{y(z)}" r="4" fill="{color}"/>')
        parts += [f'<line x1="{100+k*280}" y1="424" x2="{125+k*280}" y2="424" stroke="{color}" stroke-width="3"/>',f'<text x="{132+k*280}" y="429" font-size="12">{html.escape(label)}</text>']
    parts.append('</g></svg>');(a.out/name).write_text('\n'.join(parts))
plot('device_response.svg','35.6 µm silicon PIN detector at 850 nm','10–90% rise time (ns)',[('Simulated detector response',[float(r['device_rise_ns']) for r in rows],'#1565a5')],80)
plot('receiver_bandwidth.svg','Device speed does not set receiver bandwidth alone','−3 dB bandwidth (MHz)',[('TIA electrical (simulated)',[float(r['tia_electrical_bw_mhz']) for r in rows],'#1565a5'),('Optical receiver (calculated)',[float(r['receiver_optical_bw_mhz']) for r in rows],'#b65226')],7)
print('Wrote',a.out)
