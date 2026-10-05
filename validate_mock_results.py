"""
Check the workbooks in ./output against the ground truth saved by
generate_mock_data.py (every beat's time and type).

For each trial: true RMSSD / SDNN from the normal-to-normal intervals of the
truth (adjacent sinus beats only; the recording hole excluded) against the
pipeline's raw and corrected values, plus the Data_Quality checks that must
hold in every trial (Unique_Bad_Beats = Dropped + Corrected + Kept [D23]).
Optionally compares with workbooks from another run (e.g. the unmodified
script):  python validate_mock_results.py [other_output_dir]
"""
import os
import sys

import numpy as np
import pandas as pd

RECORDINGS = {
    'P-001_OFF': ('P-001/Equivital/EDF/P-001_OFF_truth.csv', 'P-001/Sync/SyncFile_Complete_Trial.csv'),
    'P-001_ON': ('P-001/Equivital/EDF/P-001_ON_truth.csv', 'P-001/Sync/SyncFile_Complete_Trial.csv'),
    'P-015_OFF': ('P-015/Equivital/EDF/P015_OFF_truth.csv', 'P-015/Sync/P-015_ses-01_sessionlog (2).csv'),
    'P-015_ON': ('P-015/Equivital/EDF/P015_ON_truth.csv', 'P-015/Sync/P-015_ses-01_sessionlog (2).csv'),
}
HOLES = {('P-001_OFF', 'OFF_30s_T2'): ('2026-09-15 09:11:50', '2026-09-15 09:12:00')}  # unrecorded


def truth_metrics(truth, t0, t1, hole=None):
    w = truth[(truth.clock >= t0) & (truth.clock <= t1)]
    t = w.clock.astype('int64').to_numpy() / 1e6  # ms
    k = w.kind.to_numpy()
    rr = np.diff(t)
    nn = (k[1:] == 'N') & (k[:-1] == 'N')
    if hole:
        h0, h1 = (pd.Timestamp(x).value / 1e6 for x in hole)
        nn &= ~((t[1:] > h0) & (t[:-1] < h1))
    adj = nn[1:] & nn[:-1]
    return float(np.sqrt(np.mean(np.diff(rr)[adj] ** 2))), float(np.std(rr[nn], ddof=1))


def main(out_dir='output', other=None, mock='mock_data'):
    pd.set_option('display.width', 250)
    rows, ok = [], True
    for key, (truth_f, sync_f) in RECORDINGS.items():
        wb = pd.read_excel(os.path.join(out_dir, f'HRV_Analysis_{key}.xlsx'), sheet_name=None)
        td, dq, lo = wb['Time_Domain'], wb['Data_Quality'], wb['Freq_Lomb_Scargle']
        truth = pd.read_csv(os.path.join(mock, truth_f), parse_dates=['clock'])
        truth = truth[truth.kind != 'NOISE']
        sync = pd.read_csv(os.path.join(mock, sync_f), parse_dates=['Timestamp'])
        ident = bool((dq.Unique_Bad_Beats == dq.Beats_Dropped + dq.Beats_Corrected + dq.Bad_Beats_Kept).all())
        ok &= ident
        print(f'{key}: {len(wb)} sheets, {len(td)} trials; Unique_Bad = Dropped + Corrected + Kept: {ident}')
        if other:
            old = pd.read_excel(os.path.join(other, f'HRV_Analysis_{key}.xlsx'), sheet_name='Time_Domain')
            cols = ['HR_Raw', 'SDNN_Raw', 'RMSSD_Raw', 'HR_Corrected', 'SDNN_Corrected', 'RMSSD_Corrected']
            same = np.allclose(td[cols].to_numpy(float), old[cols].to_numpy(float), equal_nan=True)
            print(f'    time-domain values identical to {other}: {same}')
        for _, r in td.iterrows():
            s = sync[sync.Trial == r.Trial]
            rm, sd = truth_metrics(truth, s.Timestamp.min(), s.Timestamp.max(), HOLES.get((key, r.Trial)))
            q = dq[dq.Trial == r.Trial].iloc[0]
            lf = lo[lo.Trial == r.Trial].iloc[0]
            rows.append({'recording': key, 'trial': r.Trial, 'beats': q.Total_Beats,
                         'bad_%': round(q.Unique_Bad_Percent, 1), 'RMSSD_true': rm, 'RMSSD_raw': r.RMSSD_Raw,
                         'RMSSD_corr': r.RMSSD_Corrected, 'SDNN_true': sd, 'SDNN_raw': r.SDNN_Raw,
                         'SDNN_corr': r.SDNN_Corrected, 'LF/HF_lomb': lf.LF_HF_Ratio_Clean,
                         'lomb_quality': lf.Quality_Clean, 'ECG': q.ECG_Channel,
                         'alternating_%': round(q.Alternating_Rhythm_Percent, 1),
                         'unrecorded_s': q.Unrecorded_Seconds})
    df = pd.DataFrame(rows)
    print()
    print(df.round(2).to_string(index=False))
    for m in ('RMSSD', 'SDNN'):
        e_raw = (df[f'{m}_raw'] / df[f'{m}_true'] - 1).abs() * 100
        e_cor = (df[f'{m}_corr'] / df[f'{m}_true'] - 1).abs() * 100
        print(f'mean |{m} error| vs truth: raw {e_raw.mean():.1f} % -> corrected {e_cor.mean():.1f} % '
              f'(median corrected {e_cor.median():.1f} %, max {e_cor.max():.1f} %)')
    short = df[df.trial.str.contains('30s')]
    print(f"30-s trials: Lomb–Scargle LF/HF empty (not 0) where not computable: {short['LF/HF_lomb'].isna().all()}")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main(other=sys.argv[1] if len(sys.argv) > 1 else None) else 1)
