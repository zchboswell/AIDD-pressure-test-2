"""Unvalidated simulated legacy converter; retain original and repair a copy."""
import math

def normalize(value, unit, relation='=', endpoint='IC50'):
    try:
        nm=float(value)*{'nM':1,'uM':1000,'M':1e9}.get(unit,1)
        p=9-math.log10(max(nm,1e-12))
        return {'status':'ok','endpoint':'IC50','value_nm':nm,
                'relation_nm':relation,'p_activity':p,'p_relation':relation}
    except (ValueError,TypeError):
        return {'status':'invalid','reason':'conversion failed'}
