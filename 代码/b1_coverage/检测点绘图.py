"""Interactive B1 bearing geometry viewer. Run: python 检测点绘图.py

Each input row: x y bearing_deg. Uses b1_coverage.py in the same folder.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import html
import json
import math
from pathlib import Path
import re
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from b1_coverage import analyze, equilateral_example

BLUE = '#245caa'
RED = '#c03f4b'
GREEN = '#158267'
INK = '#253348'
GRAY = '#8290a2'


def parse_input(count_text, coordinates, error_text):
    try:
        count = int(count_text)
    except ValueError:
        raise ValueError('点数必须为整数。') from None
    if not 1 <= count <= 100:
        raise ValueError('请输入 1～100 个检测点。')
    try:
        error = float(error_text)
    except ValueError:
        raise ValueError('误差半宽必须为数字。') from None
    if not math.isfinite(error) or not 0 < error < 90:
        raise ValueError('误差半宽必须大于 0° 且小于 90°；题面规定为 1°。')
    rows = [line.strip() for line in coordinates.splitlines() if line.strip()]
    if len(rows) != count:
        raise ValueError(f'填写的点数为 {count}，但坐标区有 {len(rows)} 行；每行应对应一个检测点。')
    observations = []
    for i, line in enumerate(rows, 1):
        parts = re.split(r'[,，;；\s]+', line)
        if len(parts) != 3:
            raise ValueError(f'第 {i} 行需要三个数：x  y  示向度。')
        try:
            x, y, bearing = map(float, parts)
        except ValueError:
            raise ValueError(f'第 {i} 行含有非数字内容。') from None
        if not all(map(math.isfinite, (x, y, bearing))):
            raise ValueError(f'第 {i} 行含有无穷大或无效数字。')
        if max(abs(x), abs(y)) > 2_000_000:
            raise ValueError(f'第 {i} 行坐标绝对值超过 2,000,000 米，请检查单位。')
        observations.append({'x': x, 'y': y, 'bearing_deg': bearing % 360})
    return observations, error


def nice_step(span, target=6):
    rough = max(span / target, 1e-12)
    power = 10 ** math.floor(math.log10(rough))
    for factor in (1, 2, 2.5, 5, 10):
        if factor * power >= rough:
            return factor * power


def fmt(value):
    return f'{value:.7g}'


def clipped_ray(s, angle, bounds):
    """Clip a forward ray to an axis-aligned rectangle in world coordinates."""
    direction = (math.cos(math.radians(angle)), math.sin(math.radians(angle)))
    lower, upper = 0.0, math.inf
    xmin, xmax, ymin, ymax = bounds
    for origin, delta, lo, hi in zip(s, direction, (xmin, ymin), (xmax, ymax)):
        if abs(delta) < 1e-14:
            if not lo <= origin <= hi:
                return None
        else:
            a, b = (lo-origin)/delta, (hi-origin)/delta
            lower, upper = max(lower, min(a, b)), min(upper, max(a, b))
    if lower > upper or not math.isfinite(upper):
        return None
    return tuple((s[0]+t*direction[0], s[1]+t*direction[1]) for t in (lower, upper))


class Scene:
    """One set of primitives for the window and an editable SVG export."""
    def __init__(self, width, height):
        self.width, self.height = width, height
        self.items = []

    def add(self, kind, coords, **options):
        self.items.append((kind, coords, options))

    def paint(self, canvas):
        canvas.delete('all')
        for kind, coords, options in self.items:
            getattr(canvas, 'create_'+kind)(*coords, **options)

    def svg(self):
        out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width}" height="{self.height}" viewBox="0 0 {self.width} {self.height}">',
               '<rect width="100%" height="100%" fill="white"/>',
               '<g font-family="Microsoft YaHei, SimHei, sans-serif">']
        for kind, coords, op in self.items:
            fill = html.escape(op.get('fill', '') or 'none', quote=True)
            stroke = html.escape(op.get('outline', ''), quote=True)
            width = op.get('width', 1)
            dash = ' stroke-dasharray="'+' '.join(map(str, op['dash']))+'"' if op.get('dash') else ''
            if kind == 'text':
                anchor = {'w':'start', 'e':'end'}.get(op.get('anchor'), 'middle')
                font = op.get('font', ('', 10))
                out.append(f'<text x="{coords[0]}" y="{coords[1]}" fill="{fill}" font-size="{font[1]*4/3}" text-anchor="{anchor}" dominant-baseline="central">{html.escape(op["text"])}</text>')
            elif kind == 'line':
                points = ' '.join(f'{coords[i]},{coords[i+1]}' for i in range(0, len(coords), 2))
                out.append(f'<polyline points="{points}" fill="none" stroke="{fill}" stroke-width="{width}"{dash}/>')
            elif kind == 'polygon':
                points = ' '.join(f'{coords[i]},{coords[i+1]}' for i in range(0, len(coords), 2))
                out.append(f'<polygon points="{points}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{width}"{dash}/>')
            elif kind == 'oval':
                x0, y0, x1, y1 = coords
                out.append(f'<ellipse cx="{(x0+x1)/2}" cy="{(y0+y1)/2}" rx="{(x1-x0)/2}" ry="{(y1-y0)/2}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{width}"{dash}/>')
            elif kind == 'rectangle':
                x0, y0, x1, y1 = coords
                out.append(f'<rect x="{x0}" y="{y0}" width="{x1-x0}" height="{y1-y0}" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{width}"/>')
        out.append('</g></svg>')
        return '\n'.join(out)


def make_scene(result, view, width, height):
    width, height = max(540, width), max(420, height)
    scene = Scene(width, height)
    scene.add('rectangle', [0, 0, width, height], fill='white', outline='white')
    bounded = result['status'] == 'bounded'
    title = '定位区域放大图' if view == 'detail' and bounded else '检测点与测向角域总览'
    scene.add('text', [width/2, 24], text=title, fill=INK, font=('Microsoft YaHei', 13, 'bold'))
    stations = [(o['x'], o['y']) for o in result['observations']]
    points = []
    if bounded:
        points.extend(result['vertices'])
        for c, radius in [(result['diameter_circle_center'], result['diameter_circle_radius_m']),
                          (result['mec_center'], result['mec_radius_m'])]:
            points.extend([(c[0]-radius, c[1]-radius), (c[0]+radius, c[1]+radius)])
    if view == 'overview' or not bounded:
        points.extend(stations)
    xmin, xmax = min(p[0] for p in points), max(p[0] for p in points)
    ymin, ymax = min(p[1] for p in points), max(p[1] for p in points)
    span = max(xmax-xmin, ymax-ymin, 1.0)
    if not bounded:
        # Visual window only: explicitly not a finite localization polygon.
        xmin -= max(100, span*.35)
        xmax += max(100, span*.35)
        ymin -= max(100, span*.35)
        ymax += max(100, span*.35)
    left, right, top, bottom = 88, width-34, 58, height-116
    scale = min((right-left)/(xmax-xmin+span*.25), (bottom-top)/(ymax-ymin+span*.25))
    cx, cy = (xmin+xmax)/2, (ymin+ymax)/2
    sx, sy = (left+right)/2, (top+bottom)/2
    def xy(p):
        return sx+(p[0]-cx)*scale, sy-(p[1]-cy)*scale
    bounds = (cx-(right-left)/(2*scale), cx+(right-left)/(2*scale),
              cy-(bottom-top)/(2*scale), cy+(bottom-top)/(2*scale))
    def line(a, b, color, line_width=1, dash=None):
        opts = {'fill': color, 'width': line_width}
        if dash:
            opts['dash'] = dash
        scene.add('line', [*xy(a), *xy(b)], **opts)
    for dim, lo, hi, target in [('x', bounds[0], bounds[1], 5), ('y', bounds[2], bounds[3], 6)]:
        step = nice_step(hi-lo, target)
        value = math.ceil(lo/step)*step
        for _ in range(30):
            if value > hi+step*1e-8:
                break
            value = 0.0 if abs(value) < step*1e-8 else value
            if dim == 'x':
                x, _ = xy((value, cy))
                scene.add('line', [x, top, x, bottom], fill='#e6ebf1')
                scene.add('text', [x, bottom+18], text=fmt(value), fill=INK, font=('Arial', 9))
            else:
                _, y = xy((cx, value))
                scene.add('line', [left, y, right, y], fill='#e6ebf1')
                scene.add('text', [left-9, y], text=fmt(value), fill=INK, anchor='e', font=('Arial', 9))
            value += step
    scene.add('rectangle', [left, top, right, bottom], outline='#9aabba', width=1)
    scene.add('text', [(left+right)/2, bottom+43], text='x / 米（东）', fill=INK, font=('Microsoft YaHei', 10))
    scene.add('text', [left, top-18], text='y / 米（北）', fill=INK, anchor='w', font=('Microsoft YaHei', 10))
    if bounded and len(result['vertices']) >= 3:
        coords = [a for p in result['vertices'] for a in xy(p)]
        scene.add('polygon', coords, fill='#e4edfc', outline=BLUE, width=2)
    for observation in result['observations']:
        s = (observation['x'], observation['y'])
        theta = observation['bearing_deg']
        for delta, dash, color in [(0, (7,4), '#7d8da5'), (-result['error_deg'], (2,4), '#aeb8c7'), (result['error_deg'], (2,4), '#aeb8c7')]:
            segment = clipped_ray(s, theta+delta, bounds)
            if segment:
                line(*segment, color, dash=dash)
    if bounded:
        for center, radius, color, dash in [(result['mec_center'], result['mec_radius_m'], GREEN, None),
                                            (result['diameter_circle_center'], result['diameter_circle_radius_m'], RED, (7,5))]:
            x, y = xy(center)
            r = radius*scale
            opts = {'outline':color, 'width':2.3}
            if dash:
                opts['dash'] = dash
            scene.add('oval', [x-r, y-r, x+r, y+r], **opts)
            scene.add('line', [x-5, y, x+5, y], fill=color, width=2)
            scene.add('line', [x, y-5, x, y+5], fill=color, width=2)
        line(*result['diameter_endpoints'], RED, 2.5)
        for i, p in enumerate(result['vertices']):
            x, y = xy(p)
            scene.add('oval', [x-3.5,y-3.5,x+3.5,y+3.5], fill=BLUE, outline=BLUE)
            if view == 'detail' and len(result['vertices']) <= 25:
                scene.add('text', [x+8,y-11], text=f'V{i+1}', fill=BLUE, anchor='w', font=('Arial', 10))
        if not result['covers'] and view == 'detail':
            x, y = xy(result['witness'])
            scene.add('oval', [x-7,y-7,x+7,y+7], outline=RED, width=2)
    if view == 'overview' or not bounded:
        for i, (s, obs) in enumerate(zip(stations, result['observations'])):
            x, y = xy(s)
            scene.add('rectangle', [x-4,y-4,x+4,y+4], fill=INK, outline=INK)
            if len(stations) <= 20:
                scene.add('text', [x+8,y-13], text=f'S{i+1} · {obs["bearing_deg"]:.3f}°', fill=INK, anchor='w', font=('Microsoft YaHei', 9))
    legends = [('定位区域', BLUE, None), ('区域直径及同直径圆', RED, (7,5)), ('最小包围圆', GREEN, None)]
    for i, (label, color, dash) in enumerate(legends):
        x = left+i*(right-left)/3
        opts = {'fill':color, 'width':2}
        if dash:
            opts['dash'] = dash
        scene.add('line', [x,height-43,x+23,height-43], **opts)
        scene.add('text', [x+29,height-43], text=label, fill=INK, anchor='w', font=('Microsoft YaHei', 9))
    if not bounded:
        note = '角域交集无界；图框仅为显示范围，不能用于计算有限直径。' if result['status']=='unbounded' else '角域交集为空：各次检测没有共同的可行位置。'
    else:
        note = f'D = {result["diameter_m"]:.6f} 米   |   R_min = {result["mec_radius_m"]:.6f} 米   |   同直径圆：'+('可覆盖' if result['covers'] else '不可覆盖')
    scene.add('text', [width/2, height-16], text=note, fill=INK, font=('Microsoft YaHei', 9))
    return scene


class Application:
    def __init__(self, root):
        self.root = root
        self.root.title('B 题问题一 · 检测点交会与圆覆盖绘图')
        self.root.geometry('1260x820')
        self.root.minsize(1050, 700)
        self.result = None
        self.job = None
        self.revision = 0
        self.executor = ThreadPoolExecutor(max_workers=1)
        self.count = tk.StringVar(value='3')
        self.error = tk.StringVar(value='1')
        self.status = tk.StringVar(value='请填写检测点，然后点击“绘图并检验”。')
        self.metrics = tk.StringVar(value='')
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)
        panel = ttk.Frame(root, padding=18)
        panel.grid(row=0,column=0,sticky='ns')
        ttk.Label(panel, text='输入检测点', font=('Microsoft YaHei', 14, 'bold')).pack(anchor='w')
        ttk.Label(panel, text='每行：x（米）  y（米）  示向度（度）', padding=(0,10,0,6)).pack(anchor='w')
        controls = ttk.Frame(panel)
        controls.pack(fill='x', pady=4)
        ttk.Label(controls, text='点数').grid(row=0,column=0,sticky='w')
        ttk.Spinbox(controls, from_=1,to=100,textvariable=self.count,width=7).grid(row=0,column=1,padx=(8,18))
        ttk.Label(controls,text='误差 ±').grid(row=0,column=2)
        ttk.Entry(controls,textvariable=self.error,width=6).grid(row=0,column=3,padx=5)
        ttk.Label(controls,text='°').grid(row=0,column=4)
        self.text = tk.Text(panel,width=37,height=19,font=('Consolas',10),wrap='none',undo=True)
        self.text.pack(fill='both',expand=True,pady=10)
        ttk.Label(panel,text='支持空格、逗号或制表符分隔。\n示向度：从正东起，逆时针计量。\n仅输入坐标不能确定定位区域。',justify='left').pack(anchor='w')
        self.button = ttk.Button(panel,text='绘图并检验',command=self.calculate)
        self.button.pack(fill='x',pady=(16,8))
        examples = ttk.Frame(panel)
        examples.pack(fill='x')
        ttk.Button(examples,text='载入三角形反例',command=lambda:self.load_example('triangle')).pack(side='left',expand=True,fill='x',padx=(0,4))
        ttk.Button(examples,text='载入可覆盖例',command=lambda:self.load_example('cover')).pack(side='left',expand=True,fill='x')
        self.export_button = ttk.Button(panel,text='导出当前图像（SVG）',command=self.export,state='disabled')
        self.export_button.pack(fill='x',pady=8)
        right = ttk.Frame(root,padding=(0,14,14,12))
        right.grid(row=0,column=1,sticky='nsew')
        right.rowconfigure(2,weight=1)
        right.columnconfigure(0,weight=1)
        ttk.Label(right,textvariable=self.status,font=('Microsoft YaHei',11,'bold'),wraplength=750).grid(row=0,column=0,sticky='ew',pady=(0,8))
        ttk.Label(right,textvariable=self.metrics,justify='left',wraplength=750).grid(row=1,column=0,sticky='ew',pady=(0,10))
        self.tabs = ttk.Notebook(right)
        self.tabs.grid(row=2,column=0,sticky='nsew')
        self.canvases = {}
        for view, label in [('detail','定位区域放大图'),('overview','检测点与测向角域总览')]:
            canvas = tk.Canvas(self.tabs,background='white',highlightthickness=0,width=740,height=620)
            self.tabs.add(canvas,text=label)
            self.canvases[view] = canvas
            canvas.bind('<Configure>',lambda event:self.draw())
        self.tabs.bind('<<NotebookTabChanged>>',lambda event:self.root.after_idle(self.draw))
        self.text.bind('<<Modified>>',self.text_changed)
        self.count.trace_add('write',self.changed)
        self.error.trace_add('write',self.changed)
        self.root.protocol('WM_DELETE_WINDOW',self.close)
        self.load_example('triangle')
        self.root.after(80,self.calculate)

    def text_changed(self, event=None):
        if self.text.edit_modified():
            self.text.edit_modified(False)
            self.changed()

    def changed(self, *_):
        self.revision += 1
        self.status.set('输入已修改；点击“绘图并检验”更新图像。')
        self.export_button.configure(state='disabled')

    def load_example(self, name):
        obs = equilateral_example() if name=='triangle' else [{'x':-1000,'y':0,'bearing_deg':0},{'x':0,'y':-1000,'bearing_deg':90}]
        self.count.set(str(len(obs)))
        self.error.set('1')
        self.text.delete('1.0','end')
        self.text.insert('1.0','\n'.join(f'{p["x"]:.14g} {p["y"]:.14g} {p["bearing_deg"]:.14g}' for p in obs))
        self.text.edit_modified(False)
        self.changed()

    def calculate(self):
        if self.job is not None:
            return
        self.text.edit_modified(False)
        try:
            observations,error = parse_input(self.count.get(),self.text.get('1.0','end'),self.error.get())
        except ValueError as exc:
            messagebox.showerror('请检查输入',str(exc),parent=self.root)
            return
        self.job_revision = self.revision
        self.status.set('正在计算定位区域、直径和最小包围圆……')
        self.button.configure(state='disabled')
        self.export_button.configure(state='disabled')
        self.job = self.executor.submit(analyze,observations,error)
        self.root.after(50,self.poll)

    def poll(self):
        if not self.job.done():
            self.root.after(50,self.poll)
            return
        job,self.job = self.job,None
        self.button.configure(state='normal')
        try:
            result = job.result()
        except Exception as exc:
            self.status.set('计算未完成，请检查输入或缩放坐标。')
            messagebox.showerror('计算错误',str(exc),parent=self.root)
            return
        if self.job_revision != self.revision:
            self.status.set('计算期间输入已修改，请重新点击“绘图并检验”。')
            return
        self.show_result(result)

    def show_result(self,result):
        self.result = result
        if result['status']=='bounded':
            self.status.set('同直径圆可以覆盖整个定位区域。' if result['covers'] else '同直径圆不能覆盖定位区域，调整圆心也无法实现覆盖。')
            self.metrics.set(f'区域顶点：{result["vertex_count"]} 个    区域直径 D：{result["diameter_m"]:.6f} 米\n'
                             f'同直径圆半径 D/2：{result["diameter_circle_radius_m"]:.6f} 米    最小包围圆半径：{result["mec_radius_m"]:.6f} 米\n'
                             f'最大超出距离：{max(0,result["maximum_radial_excess_m"]):.6f} 米    判定距离容差：{result["tolerance_m"]:g} 米')
        else:
            self.status.set('角域交集无界：这些检测点不能确定有限的定位多边形。' if result['status']=='unbounded' else '角域交集为空：请检查坐标、示向度及误差。')
            self.metrics.set('没有可计算的有限区域直径；显示测向射线供检查。图框仅用于显示。')
            self.tabs.select(1)
        self.export_button.configure(state='normal')
        self.draw()

    def current_scene(self):
        view = 'detail' if self.tabs.index('current')==0 else 'overview'
        canvas = self.canvases[view]
        return make_scene(self.result,view,canvas.winfo_width(),canvas.winfo_height())

    def draw(self):
        if self.result:
            view = 'detail' if self.tabs.index('current')==0 else 'overview'
            self.current_scene().paint(self.canvases[view])

    def export(self):
        if not self.result:
            return
        path = filedialog.asksaveasfilename(parent=self.root,title='保存图像',defaultextension='.svg',initialfile='B1_定位区域.svg',filetypes=[('SVG 矢量图','*.svg')])
        if path:
            try:
                Path(path).write_text(self.current_scene().svg(),encoding='utf-8')
            except OSError as exc:
                messagebox.showerror('保存失败',str(exc),parent=self.root)
                return
            messagebox.showinfo('已保存','图像已保存，可用浏览器打开或插入论文。',parent=self.root)

    def close(self):
        self.executor.shutdown(wait=False,cancel_futures=True)
        self.root.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,help='与 b1_coverage.py 相同的 JSON 输入')
    parser.add_argument('--svg',type=Path,help='无需启动窗口，直接输出 SVG 图像')
    parser.add_argument('--view',choices=['detail','overview'],default='detail')
    args = parser.parse_args()
    if args.svg:
        data = json.loads(args.input.read_text(encoding='utf-8-sig')) if args.input else {'observations':equilateral_example()}
        # Apply the same validation to the noninteractive entry point.
        obs,error = parse_input(str(len(data['observations'])), '\n'.join(f'{o["x"]} {o["y"]} {o["bearing_deg"]}' for o in data['observations']),str(data.get('error_deg',1)))
        result = analyze(obs,error)
        args.svg.parent.mkdir(parents=True,exist_ok=True)
        args.svg.write_text(make_scene(result,args.view,1100,800).svg(),encoding='utf-8')
        print(str(args.svg))
    else:
        root = tk.Tk()
        Application(root)
        root.mainloop()


if __name__=='__main__':
    main()
