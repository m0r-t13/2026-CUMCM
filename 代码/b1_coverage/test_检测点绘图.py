"""Input, export and real Tk canvas smoke checks for the drawing program."""
import math
from pathlib import Path
import tempfile
import time
import tkinter as tk
import unittest
import xml.etree.ElementTree as ET

from b1_coverage import analyze, equilateral_example
from 检测点绘图 import Application, clipped_ray, make_scene, parse_input


class ViewerTests(unittest.TestCase):
    def test_input(self):
        obs, err = parse_input('2', '-1000，0，360\n0\t-1000\t90', '1')
        self.assertEqual(obs[0]['bearing_deg'], 0)
        self.assertEqual(err, 1)
        for count, data, error in [('2','0 0 1','1'), ('1','0 0','1'),
                                    ('1','nan 0 90','1'), ('1','0 0 90','90'),
                                    ('1','0 0 90','nan'), ('0','','1')]:
            with self.assertRaises(ValueError):
                parse_input(count,data,error)

    def test_rays(self):
        self.assertEqual(clipped_ray((-2,0),0,(-1,1,-1,1)),((-1,0),(1,0)))
        self.assertIsNone(clipped_ray((2,0),0,(-1,1,-1,1)))

    def test_svg_and_equal_scale(self):
        r = analyze(equilateral_example())
        for view in ['detail','overview']:
            scene = make_scene(r,view,1100,800)
            doc = ET.fromstring(scene.svg())
            ellipses = doc.findall('.//{http://www.w3.org/2000/svg}ellipse')
            self.assertGreater(len(ellipses),2)
            for circle in ellipses:
                self.assertAlmostEqual(float(circle.attrib['rx']),float(circle.attrib['ry']),places=8)
            for _,coords,_ in scene.items:
                self.assertTrue(all(math.isfinite(x) for x in coords))

    def test_window_calculation_and_redraw(self):
        root = tk.Tk()
        root.withdraw()
        app = Application(root)
        deadline = time.monotonic()+15
        while app.result is None and time.monotonic()<deadline:
            root.update()
            time.sleep(.02)
        self.assertIsNotNone(app.result)
        self.assertFalse(app.result['covers'])
        self.assertIn('不能覆盖',app.status.get())
        for view,canvas in app.canvases.items():
            make_scene(app.result,view,900,700).paint(canvas)
            self.assertGreater(len(canvas.find_all()),20)
        app.load_example('cover')
        app.calculate()
        while app.job is not None and time.monotonic()<deadline:
            root.update()
            time.sleep(.02)
        self.assertIsNone(app.job)
        self.assertTrue(app.result['covers'])
        for obs,status in [([{'x':0,'y':0,'bearing_deg':0}],'unbounded'),
                           ([{'x':0,'y':0,'bearing_deg':180},{'x':100,'y':0,'bearing_deg':0}],'empty')]:
            result = analyze(obs)
            self.assertEqual(result['status'],status)
            app.show_result(result)
            root.update()
            self.assertGreater(len(app.canvases['overview'].find_all()),10)
        app.close()


if __name__=='__main__':
    unittest.main(verbosity=2)
