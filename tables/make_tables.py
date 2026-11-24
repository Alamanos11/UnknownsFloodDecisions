"""Paper 1 – builds the manuscript tables as CSV (T1, T3, T4; T2 from T2_who_decides.md; T5 by triage_costs.py).
Numbers from the case READMEs / METHODS.md (T1, T3) and from the raw result CSVs in raw/ (T4). usage: python make_tables.py"""
import os, csv, re
H = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(H, 'raw')
def write(name, rows):
    with open(os.path.join(H, name), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print('written', name)
# T1 – torrents and models (sources: case READMEs, 03_hecras READMEs, serres_programme_status)
T1 = [
 dict(case='SER-01', stream='Agioi Anargyroi', catchment_km2=76.29, setting='urban fan; channel on a ridge through Serres', outlet_E=462038, outlet_N=4550592,
      Q_T50_m3s=149.2, Q_T100_m3s=204.1, Q_T1000_m3s=470.8, domain_km2=9.92, cells=102206, C1_T50_flooded_ha=227, buildings_in_domain=5882,
      flood_records='18/06/2004; 21/06/2020 (73 buildings, FRMP EL11)'),
 dict(case='SER-02', stream='Agia Varvara (Kamenitsa)', catchment_km2=32.71, setting='urban edge; regulated reach (PKM) to the plain', outlet_E=460208, outlet_N=4548692,
      Q_T50_m3s=67.0, Q_T100_m3s=91.0, Q_T1000_m3s=207.8, domain_km2=11.78, cells=122126, C1_T50_flooded_ha=137, buildings_in_domain=445,
      flood_records='21/06/2020 (overflow, no location; FRMP EL11)'),
 dict(case='SER-03', stream='Agios Ioannis', catchment_km2=8.77, setting='spring-fed; sheet flow on a fan; EO12 embankment', outlet_E=465118, outlet_N=4547957,
      Q_T50_m3s=20.0, Q_T100_m3s=26.3, Q_T1000_m3s=56.0, domain_km2=6.14, cells=60622, C1_T50_flooded_ha=67.8, buildings_in_domain=824,
      flood_records='none recorded (1925 only)'),
 dict(case='SER-04', stream='Krousovitis', catchment_km2=270.42, setting='gorge -> Sidirokastro -> plain -> Belitsa canal', outlet_E=448748, outlet_N=4565072,
      Q_T50_m3s=429.0, Q_T100_m3s=579.2, Q_T1000_m3s=1306.5, domain_km2=55.35, cells=101517, C1_T50_flooded_ha=856, buildings_in_domain=2884,
      flood_records='11/01/2021 (FRMP EL11); 26/01/2026 (state of emergency)'),
]
write('T1_torrents_models.csv', T1)
# T2 – from markdown table
md = open(os.path.join(H, 'T2_who_decides.md'), encoding='utf-8').read().splitlines()
rows = [l for l in md if l.startswith('| **D')]
T2 = []
for l in rows:
    c = [x.strip() for x in l.strip('|').split('|')]
    T2.append(dict(decision=re.sub(r'\*', '', c[0]), who_decides_legal_basis=c[1], evidence_today=c[2], our_metric=c[3]))
write('T2_who_decides.csv', T2)
# T3 – factorial of unknowns and levers
T3 = [
 dict(factor='channel geometry', type='epistemic unknown', levels='C0 DEM only (crossings opened); C1 assumed channel; C+ larger channel', runs_per_torrent='all designs'),
 dict(factor='antecedent soil moisture', type='epistemic unknown', levels='AMC II (CN II); AMC III (CN III)', runs_per_torrent='3 (T50 CN III on C0/C1/C+)'),
 dict(factor='design return period', type='regulatory choice', levels='T50; T100; T1000 (C1 only)', runs_per_torrent='3 + 1'),
 dict(factor='maintenance state', type='management lever', levels='base n 0.04; overgrown n 0.08; cleared n 0.03; half-blocked crossings; bed +0.5 m', runs_per_torrent='8 (4 states x C1/C+)'),
 dict(factor='inflow ramp (D3)', type='threshold design', levels='15 steps 5-100 % of Q1000, hold 60-120 min', runs_per_torrent='3 (C0/C1/C+)'),
 dict(factor='events', type='plausibility', levels='SER-01 2004, 2020; SER-04 Jan 2021 A (gauge, CERRA), B', runs_per_torrent='SER-01 2; SER-04 3'),
 dict(factor='total HEC-RAS 2D runs', type='', levels='46 core and event runs + 44 maintenance and stepped-inflow runs = 90', runs_per_torrent=''),
]
write('T3_factorial.csv', T3)
# T4 – decision answers (synthesis of the results; numbers from raw/)
T4 = [
 dict(decision='D1 which stream / reach first', decidable_without_survey='yes', what_changes_it='little: stream order fixed in 17/18 uniform scenarios (18/18 for buildings >= 0.50 m); SER-01 first in 93-100 % of mixed states; reach rho 0.64-0.98',
      cheapest_information='none needed to rank streams; reach ranking in SER-04 unstable'),
 dict(decision='D2 is upkeep worth it, where', decidable_without_survey='yes, per torrent', what_changes_it='maintenance state as large as the channel unknown (line position: maintenance 39-52 % vs channel 21-33 % in SER-01/02/04)',
      cheapest_information='model each cleaning before contracting; SER-03 cleaning has no effect (-2...+6 %)'),
 dict(decision='D3 at what discharge damage starts', decidable_without_survey='as a band', what_changes_it='channel shifts Q for 50 buildings by a factor of about 2 between adjacent states (up to ~4 in SER-02; upper bounds in SER-04); soil moisture shifts the warning rain by 35-45 %',
      cheapest_information='channel survey at the controlling reach + soil-moisture indicator'),
 dict(decision='D4 where the legal line falls', decidable_without_survey='only SER-03, after a check of the culverts', what_changes_it='channel (69-75 % of line position with C0; 37-52 % between C1 and C+); T50 vs T100 11-20 % (27-48 % between C1 and C+); band resolvable by a survey 44-1 548 ha with 29-3 174 buildings',
      cheapest_information='targeted channel survey: EUR 13-20 per building at stake in SER-01, EUR 120-195 in SER-04'),
]
write('T4_decision_answers.csv', T4)
