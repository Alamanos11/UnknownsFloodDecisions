"""Paper 1, Table 5 – indicative cost of the information that would settle D4, per case, set against the buildings whose
in/out status depends on the unknowns (D4 band). Indicative only (regulatory fee scale, 2025 prices, no VAT, no contingencies).

Inputs
  D4_band_summary.csv (02_CASES/SER-P1_paper_decisions/D4_floodlines): buildings 'sometimes' flooded over 14 plausible T50 states
  modelled channel length per case (centrelines used to build the models; measured 23/11/2026)
Cost items (sources in 04_PAPERS/shared/costs_sources.md)
  targeted channel survey: 1:500 ground survey of a 100 m corridor along the modelled channel (= 100 stremmata per km),
     unbuilt EUR 30 / built EUR 60 per stremma (TOP.5/1.1, TOP.6/2.2) x tk 1.435  -> EUR 4.3-8.6 k per km
  hydraulic delineation fee (YDR.14, beta = 3 full delineation): A = 60*beta*(5 + 20 L^(2/3) + 2.5 F^(1/3)) * tk
  full study budgets of real tenders: EUR 0.12-0.75 M (Arta, Spetses, Nisoi, Chiliadou) - range only
usage: python triage_costs.py <D4_band_summary.csv> <out.csv>
"""
import sys, csv
TK = 1.435
CASE = {  # modelled channel length (km), catchment area (km2)
    'SER-01': (5.04, 76.29), 'SER-02': (6.69, 32.71), 'SER-03': (3.10, 8.77), 'SER-04': (13.18, 270.42)}
band = {r['case']: r for r in csv.DictReader(open(sys.argv[1]))}
rows = []
for c, (L, F) in CASE.items():
    surv_lo, surv_hi = L * 100 * 30 * TK, L * 100 * 60 * TK
    ydr14 = 60 * 3 * (5 + 20 * L ** (2 / 3) + 2.5 * F ** (1 / 3)) * TK
    b = int(band[c]['bld_sometimes']); a = float(band[c]['A_sometimes_ha'])
    rows.append(dict(case=c, channel_km=L, catchment_km2=F, band_ha=round(a), band_buildings=b,
                     buildings_per_km=round(b / L, 1),
                     survey_EUR_lo=round(surv_lo, -2), survey_EUR_hi=round(surv_hi, -2),
                     hydraulic_fee_YDR14_EUR=round(ydr14, -2),
                     EUR_per_building_lo=round((surv_lo + ydr14) / b) if b else '',
                     EUR_per_building_hi=round((surv_hi + ydr14) / b) if b else '',
                     full_study_EUR_per_building_lo=round(120000 / b) if b else '',
                     full_study_EUR_per_building_hi=round(750000 / b) if b else ''))
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
