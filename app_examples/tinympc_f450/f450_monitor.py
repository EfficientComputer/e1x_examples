#!/usr/bin/env python3
"""
f450_monitor.py — live web dashboard for the F450 sensor-check firmware
(tinympc_f450_check). Reads the CSV telemetry over serial, serves a real-time
3D attitude view + calibration readouts at http://localhost:8420.

  python3 f450_monitor.py -p /dev/cu.usbmodem105 -b 108000   # live
  python3 f450_monitor.py --sim                              # synthetic (no HW)

CSV: F450CHK,step,gyro_mrads[3],anorm_mg[3],mag[3],baro_mm,gps_fix,gpsx_mm,
     gpsy_mm,roll_md,pitch_md,yaw_md,px_mm,py_mm,pz_mm,vx,vy,vz,bgx,bgy,bgz
"""
import argparse, json, math, threading, time, http.server, socketserver

FIELDS = ["step","gx","gy","gz","ax","ay","az","mx","my","mz","baro",
          "gfix","gpx","gpy","roll","pitch","yaw","px","py","pz","vx","vy","vz","bgx","bgy","bgz"]
latest = {k: 0 for k in FIELDS}
latest["conn"] = 0

def parse_line(line):
    p = line.strip().split(",")
    if len(p) < 27 or p[0] != "F450CHK":
        return None
    try:
        vals = [int(x) for x in p[1:27]]
    except ValueError:
        return None
    return dict(zip(FIELDS, vals))

def serial_reader(port, baud):
    import serial
    while True:
        try:
            s = serial.Serial(port, baud, timeout=1); s.setDTR(False); s.setRTS(False)
            latest["conn"] = 1
            for raw in s:
                d = parse_line(raw.decode("ascii", "ignore"))
                if d: latest.update(d)
        except Exception as e:
            latest["conn"] = 0; print("serial:", e); time.sleep(1)

def sim_reader():
    """synthetic telemetry: sweep attitude + a vertical bounce so the GUI can be
    verified with no hardware."""
    latest["conn"] = 2; t = 0.0
    while True:
        t += 0.05
        latest.update(step=int(t*20),
            roll=int(20000*math.sin(t*0.7)), pitch=int(12000*math.sin(t*0.5)),
            yaw=int(((t*20) % 360)*1000),
            gx=int(300*math.cos(t*0.7)), gy=int(200*math.cos(t*0.5)), gz=200,
            ax=int(30*math.sin(t*0.7)), ay=int(-30*math.sin(t*0.5)), az=1000,
            mx=1200, my=int(400*math.sin(t*0.3)), mz=-2100,
            baro=int(500+300*math.sin(t*0.4)),
            gfix=1, gpx=int(200*math.sin(t*0.2)), gpy=int(150*math.cos(t*0.2)),
            px=int(180*math.sin(t*0.2)), py=int(140*math.cos(t*0.2)), pz=int(500+300*math.sin(t*0.4)),
            vx=int(50*math.cos(t*0.2)), vy=0, vz=int(40*math.cos(t*0.4)),
            bgx=2, bgy=-3, bgz=1)
        time.sleep(0.05)

PAGE = """<!doctype html><html><head><meta charset=utf-8><title>F450 monitor</title>
<style>
 body{margin:0;background:#0b0e14;color:#c8d3e0;font:13px/1.4 "SF Mono",Menlo,monospace}
 #wrap{display:flex;height:100vh}#view{flex:1;position:relative}
 #hud{width:320px;padding:14px 16px;overflow:auto;border-left:1px solid #1e2836}
 h1{font-size:14px;color:#5ec8ff;margin:0 0 10px}h2{font-size:11px;color:#7a8aa0;margin:14px 0 4px;text-transform:uppercase;letter-spacing:.08em}
 .big{font-size:26px;color:#fff}.row{display:flex;justify-content:space-between;padding:2px 0}
 .k{color:#7a8aa0}.v{color:#e6edf5}.warn{color:#ffb454}.ok{color:#7ee787}
 .bar{height:6px;background:#1e2836;border-radius:3px;overflow:hidden;margin:3px 0}
 .bar>i{display:block;height:100%;background:#5ec8ff}
 #conn{position:absolute;top:10px;left:12px;padding:3px 8px;border-radius:4px;background:#1e2836}
</style></head><body><div id=wrap>
<div id=view><canvas id=c></canvas><div id=conn>…</div></div>
<div id=hud>
 <h1>F450 SENSOR CHECK</h1>
 <h2>Attitude</h2>
 <div class=row><span class=k>roll</span><span class="v big" id=roll>0</span></div>
 <div class=row><span class=k>pitch</span><span class="v big" id=pitch>0</span></div>
 <div class=row><span class=k>yaw</span><span class="v big" id=yaw>0</span></div>
 <h2>Position (m) / Altitude</h2>
 <div class=row><span class=k>x n</span><span class=v id=px>0</span></div>
 <div class=row><span class=k>y e</span><span class=v id=py>0</span></div>
 <div class=row><span class=k>z up</span><span class=v id=pz>0</span></div>
 <div class=row><span class=k>vel</span><span class=v id=vel>0</span></div>
 <h2>Sensors</h2>
 <div class=row><span class=k>gyro mrad/s</span><span class=v id=gyro>0</span></div>
 <div class=row><span class=k>accel /g</span><span class=v id=acc>0</span></div>
 <div class=row><span class=k>|accel| (≈1 level)</span><span class=v id=amag>0</span></div>
 <div class=row><span class=k>mag</span><span class=v id=mag>0</span></div>
 <div class=row><span class=k>baro alt</span><span class=v id=baro>0</span></div>
 <div class=row><span class=k>GPS</span><span class=v id=gps>—</span></div>
 <div class=row><span class=k>gyro bias</span><span class=v id=bg>0</span></div>
</div></div>
<script type="importmap">{"imports":{"three":"https://cdn.jsdelivr.net/npm/three@0.163.0/build/three.module.js"}}</script>
<script type=module>
import * as THREE from 'three';
const cv=document.getElementById('c');
const rn=new THREE.WebGLRenderer({canvas:cv,antialias:true});
const sc=new THREE.Scene(); sc.background=new THREE.Color(0x0b0e14);
const vw=()=>document.getElementById('view');
const cam=new THREE.PerspectiveCamera(50,1,0.1,100); cam.position.set(0.6,0.5,0.8);
cam.lookAt(0,0,0);
sc.add(new THREE.AmbientLight(0xffffff,0.7));
const dl=new THREE.DirectionalLight(0xffffff,0.8); dl.position.set(1,2,1); sc.add(dl);
sc.add(new THREE.GridHelper(2,10,0x223,0x162030));
// quad: body + 4 arms + rotors (X config)
const quad=new THREE.Group();
const body=new THREE.Mesh(new THREE.BoxGeometry(0.12,0.04,0.12),new THREE.MeshStandardMaterial({color:0x30506f}));
quad.add(body);
const armMat=new THREE.MeshStandardMaterial({color:0x9aa7b4});
const rotMat=[0x5ec8ff,0x5ec8ff,0xffb454,0xffb454];
const pos=[[0.18,0,0.18],[-0.18,0,-0.18],[0.18,0,-0.18],[-0.18,0,0.18]];
pos.forEach((p,i)=>{
  const arm=new THREE.Mesh(new THREE.CylinderGeometry(0.008,0.008,0.26),armMat);
  arm.position.set(p[0]/2,0,p[2]/2); arm.rotation.z=Math.PI/2;
  arm.lookAt(new THREE.Vector3(p[0],0,p[2])); quad.add(arm);
  const r=new THREE.Mesh(new THREE.CylinderGeometry(0.07,0.07,0.006,20),new THREE.MeshStandardMaterial({color:rotMat[i],transparent:true,opacity:0.55}));
  r.position.set(p[0],0.02,p[2]); quad.add(r);
});
// nose marker (+x forward)
const nose=new THREE.Mesh(new THREE.ConeGeometry(0.02,0.06,12),new THREE.MeshStandardMaterial({color:0xff5e7e}));
nose.position.set(0.09,0,0); nose.rotation.z=-Math.PI/2; quad.add(nose);
sc.add(quad);
function resize(){const w=vw().clientWidth,h=vw().clientHeight;rn.setSize(w,h,false);cam.aspect=w/h;cam.updateProjectionMatrix();}
new ResizeObserver(resize).observe(vw()); resize();
let att={roll:0,pitch:0,yaw:0};
function loop(){
  // roll about +x(fwd), pitch about +z(right→ our z), yaw about +y(up)
  quad.rotation.set(0,0,0);
  quad.rotateY(att.yaw); quad.rotateX(att.pitch); quad.rotateZ(att.roll);
  rn.render(sc,cam); requestAnimationFrame(loop);
} loop();
const $=id=>document.getElementById(id);
const es=new EventSource('/stream');
es.onmessage=e=>{
  const d=JSON.parse(e.data);
  const r=d.roll/1000,p=d.pitch/1000,y=d.yaw/1000;
  att={roll:r*Math.PI/180,pitch:p*Math.PI/180,yaw:-y*Math.PI/180};
  $('roll').textContent=r.toFixed(1)+'°';$('pitch').textContent=p.toFixed(1)+'°';$('yaw').textContent=y.toFixed(1)+'°';
  $('px').textContent=(d.px/1000).toFixed(2);$('py').textContent=(d.py/1000).toFixed(2);$('pz').textContent=(d.pz/1000).toFixed(2);
  $('vel').textContent=`${(d.vx/1000).toFixed(2)}, ${(d.vy/1000).toFixed(2)}, ${(d.vz/1000).toFixed(2)} m/s`;
  $('gyro').textContent=`${d.gx}, ${d.gy}, ${d.gz}`;
  const am=Math.sqrt(d.ax*d.ax+d.ay*d.ay+d.az*d.az)/1000;
  $('acc').textContent=`${d.ax}, ${d.ay}, ${d.az} mg`;
  const amEl=$('amag'); amEl.textContent=am.toFixed(3)+' g'; amEl.className='v '+(Math.abs(am-1)<0.1?'ok':'warn');
  $('mag').textContent=`${d.mx}, ${d.my}, ${d.mz}`;
  $('baro').textContent=(d.baro/1000).toFixed(2)+' m';
  const g=$('gps'); g.textContent=d.gfix?`3D  N ${(d.gpx/1000).toFixed(1)} E ${(d.gpy/1000).toFixed(1)} m`:'no fix'; g.className='v '+(d.gfix?'ok':'warn');
  $('bg').textContent=`${d.bgx}, ${d.bgy}, ${d.bgz} mrad/s`;
  const c=document.getElementById('conn');
  c.textContent=d.conn==2?'SIM':(d.conn?'LIVE':'no link'); c.style.color=d.conn?'#7ee787':'#ffb454';
};
</script></body></html>"""

class H(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def do_GET(self):
        if self.path == "/":
            b = PAGE.encode(); self.send_response(200)
            self.send_header("Content-Type","text/html"); self.send_header("Content-Length",str(len(b)))
            self.end_headers(); self.wfile.write(b)
        elif self.path == "/stream":
            self.send_response(200); self.send_header("Content-Type","text/event-stream")
            self.send_header("Cache-Control","no-cache"); self.end_headers()
            try:
                while True:
                    self.wfile.write(b"data: "+json.dumps(latest).encode()+b"\n\n"); self.wfile.flush()
                    time.sleep(0.05)
            except Exception: pass
        else:
            self.send_error(404)

class Srv(socketserver.ThreadingMixIn, http.server.HTTPServer): daemon_threads=True

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-p","--port", default="/dev/cu.usbmodem105")
    ap.add_argument("-b","--baud", type=int, default=108000)
    ap.add_argument("--http", type=int, default=8420)
    ap.add_argument("--sim", action="store_true", help="synthetic telemetry, no hardware")
    a = ap.parse_args()
    threading.Thread(target=sim_reader if a.sim else serial_reader,
                     args=() if a.sim else (a.port, a.baud), daemon=True).start()
    print(f"F450 monitor → http://localhost:{a.http}   ({'SIM' if a.sim else a.port})")
    Srv(("127.0.0.1", a.http), H).serve_forever()
