from module1 import Design
foil_thickness = 0.12
def mass_builtup(b,mac):
    s = mac*b
    m_LE = 0.15*s*0.58*2
    m_spar = ((0.0254*0.25)*(b*mac*foil_thickness)*100)+(2*b*0.0254*0.25*0.002*1750)
    m_ribs = 0.85*s*0.5*foil_thickness*(100*0.025)
    m_TE = b*1e-6*200
    m_skin = s*0.75*2*0.036584
    return m_LE,m_spar,m_ribs,m_TE,m_skin,sum([m_LE,m_spar,m_ribs,m_TE,m_skin])

print(mass_builtup(2.62,0.294))