import sys
sys.path.insert(0, r"C:/program1/Program/2019.2vivadoprj-master/DSB/analysis")
import design_n12_unit
pts=design_n12_unit.layout_rect_3x4(20,20)
cases=[(0,0),(15,0)]
max_sl=-1e9
p99_sl=-1e9
axis, tx, ty, theta_deg, dirs, elem, mask = design_n12_unit.GRID
for sx,sy in cases:
    w=design_n12_unit.steering_vector(pts,sx,sy)
    pts_m = design_n12_unit.np.column_stack([pts[:,0]*1e-3, pts[:,1]*1e-3, design_n12_unit.np.zeros(len(pts))])
    phase = design_n12_unit.np.exp(1j*design_n12_unit.K*design_n12_unit.np.einsum('ijk,kl->ijl', dirs, pts_m.T))
    resp = elem*design_n12_unit.np.dot(phase, w)
    target_amp = abs(design_n12_unit.np.dot(design_n12_unit.np.exp(1j*design_n12_unit.K*design_n12_unit.np.dot(pts_m, design_n12_unit.np.array([design_n12_unit.math.sin(design_n12_unit.math.radians(sx)), design_n12_unit.math.sin(design_n12_unit.math.radians(sy)), design_n12_unit.math.sqrt(max(0.0,1.0-design_n12_unit.math.sin(design_n12_unit.math.radians(sx))**2-design_n12_unit.math.sin(design_n12_unit.math.radians(sy))**2))]))), w))
    db = 20.0*design_n12_unit.np.log10(design_n12_unit.np.abs(resp)/max(target_amp,1e-18)+1e-12)
    db[~mask] = design_n12_unit.np.nan
    sep=design_n12_unit.angular_separation_deg(tx,ty,sx,sy)
    outside=(sep>11.0)&mask
    sl=float(design_n12_unit.np.nanmax(db[outside]))
    p99=float(design_n12_unit.np.nanpercentile(db[outside],99.0))
    max_sl=max(max_sl,sl)
    p99_sl=max(p99_sl,p99)
print(max_sl,p99_sl)
