"""
Synthetic ECG test data for the HRV pipeline (`main`).

Writes EDF recordings, sync logs and trial run sheets under ./mock_data, laid
out exactly as the pipeline's PARTICIPANTS registry expects
(<BASE>/<PID>/Equivital/EDF/..., <BASE>/<PID>/Sync/...), for two registry
participants whose entries exercise the special cases of the registry:

  P-001  OFF: two EDF files with a 10-s recording hole inside OFF_30s_T2 [D26].
         ON : the 5-min sitting is in a separate 'Rest' file, matched by the
              second pattern of the ON entry [D25].
  P-015  OFF: ECG 1 is noisy (electrode artefact, dropouts) during
              OFF_5min_Sitting, ECG 2 is clean - the registry reads that trial
              from ECG 2 [D27]. The OFF recording also covers an 'ON_' trial
              that must be skipped [D20].
         ON : a sustained 60-s bigeminy run in ON_5min_Standing [D21/D28] and a
              2-s placeholder trial in the sync log.

Each recording contains a physiological sinus rhythm (respiratory sinus
arrhythmia ~0.25 Hz, Mayer waves ~0.1 Hz, random beat-to-beat variability),
baseline wander (respiration + slow drift), mains interference and broadband
noise, plus injected events with known times:
  PVC   premature ventricular beat (wide QRS, no P, inverted T), coupling
        ~65 % of the RR, full compensatory pause
  APC   premature atrial beat (normal QRS), coupling ~75 %, sinus reset
  MISS  a beat with a tiny R wave (dropped by the detector -> missed beat)
  NOISE 0.6-s burst of QRS-band noise (spurious detections -> extra beats)
The truth (every beat time and type) is saved next to the EDFs as
*_truth.csv so the reconstruction can be checked.

Only numpy / scipy / mne (+ edfio, MNE's EDF writer) are used.
Run:  python generate_mock_data.py            (seeded, reproducible)
"""
import os
from datetime import datetime, timedelta, timezone

import mne
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mock_data')
FS = 256.0
SESSION_DATE = datetime(2026, 9, 15, tzinfo=timezone.utc)
RNG = np.random.default_rng(20261005)


# ─── waveform ───────────────────────────────────────────────────────────────
def _gauss(t, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((t - mu) / sigma) ** 2)


def beat_wave(t, kind, amp_scale=1.0):
    """One PQRST complex (mV) on t (s, relative to the R peak)."""
    if kind == 'PVC':  # ventricular: no P wave, wide bizarre QRS, discordant T
        return (_gauss(t, 1.35, 0.0, 0.030) + _gauss(t, -0.55, 0.06, 0.030)
                + _gauss(t, -0.40, 0.30, 0.070)) * amp_scale
    w = (_gauss(t, 0.12, -0.20, 0.025) + _gauss(t, -0.10, -0.030, 0.010)
         + _gauss(t, 1.00, 0.0, 0.011) + _gauss(t, -0.25, 0.032, 0.012)
         + _gauss(t, 0.30, 0.26, 0.050))
    if kind == 'APC':  # premature atrial: different P wave, normal QRS
        w = w - _gauss(t, 0.12, -0.20, 0.025) + _gauss(t, -0.06, -0.16, 0.02)
    return w * amp_scale


# ─── rhythm ─────────────────────────────────────────────────────────────────
def sinus_beats(t0, t1, mean_rr, rsa=0.035, mayer=0.020, jitter=0.012):
    """Sinus beat times in [t0, t1): RR modulated by RSA (0.25 Hz) and Mayer
    waves (0.1 Hz) plus random variability (all in s)."""
    ph_r, ph_m = RNG.uniform(0, 2 * np.pi, 2)
    times, t = [], t0 + RNG.uniform(0, mean_rr)
    while t < t1:
        times.append(t)
        rr = (mean_rr + rsa * np.sin(2 * np.pi * 0.25 * t + ph_r)
              + mayer * np.sin(2 * np.pi * 0.10 * t + ph_m) + RNG.normal(0, jitter))
        t += max(rr, 0.35)
    return times


def build_beats(blocks, events):
    """blocks: [(t0, t1, mean_rr)] of sinus rhythm. events: [(time, kind)]
    applied to the sinus beat nearest after `time`.
    Returns a DataFrame of beats: time, kind, amp."""
    beats = []
    for t0, t1, rr in blocks:
        beats += [[t, 'N', 1.0] for t in sinus_beats(t0, t1, rr)]
    beats.sort(key=lambda b: b[0])
    times = np.array([b[0] for b in beats])
    for t_ev, kind in events:
        i = int(np.searchsorted(times, t_ev))
        if i < 2 or i >= len(beats) - 2:
            continue
        prev_rr = beats[i - 1][0] - beats[i - 2][0]
        if kind == 'PVC':  # early beat, next sinus beat on time (compensatory)
            beats[i][0] = beats[i - 1][0] + 0.65 * prev_rr
            beats[i][1] = 'PVC'
            beats[i + 1][0] = beats[i - 1][0] + 2.0 * prev_rr
        elif kind == 'APC':  # early beat that resets the sinus node
            new_t = beats[i - 1][0] + 0.75 * prev_rr
            shift = new_t - beats[i][0]
            for b in beats[i:]:
                b[0] += shift
            beats[i][1] = 'APC'
        elif kind == 'MISS':  # beat present, R wave too small to detect
            beats[i][1], beats[i][2] = 'MISS', 0.08
        elif kind == 'NOISE':
            beats.append([t_ev, 'NOISE', 0.0])
        times = np.array([b[0] for b in beats])
    return pd.DataFrame(beats, columns=['time_s', 'kind', 'amp']).sort_values('time_s').reset_index(drop=True)


def bigeminy(beats, t0, t1, coupling=0.80):
    """Every other beat in [t0, t1) becomes an early (atrial) beat with a
    resetting pause, giving a sustained short-long alternation (as in P-014)."""
    b = beats.copy()
    idx = b.index[(b.time_s >= t0) & (b.time_s < t1) & (b.kind == 'N')].to_numpy()
    for j, i in enumerate(idx[1:-1:2]):
        prev = b.at[i - 1, 'time_s']
        rr = b.at[i, 'time_s'] - prev
        b.at[i, 'time_s'] = prev + coupling * rr
        b.at[i, 'kind'] = 'APC'
    return b.sort_values('time_s').reset_index(drop=True)


# ─── signal ─────────────────────────────────────────────────────────────────
def synth_ecg(beats, duration, lead_gain=1.0, noise_sd=0.02):
    n = int(round(duration * FS))
    t = np.arange(n) / FS
    ecg = np.zeros(n)
    half = int(0.45 * FS)
    for bt, kind, amp in beats[['time_s', 'kind', 'amp']].itertuples(index=False):
        if kind == 'NOISE':  # QRS-band burst (spurious detections)
            i0 = int(bt * FS)
            seg = RNG.normal(0, 0.6, int(0.6 * FS))
            seg = sosfiltfilt(butter(2, [8, 25], 'bandpass', fs=FS, output='sos'), seg) * 3
            ecg[i0:i0 + len(seg)] += seg[:max(0, min(len(seg), n - i0))]
            continue
        c = int(round(bt * FS))
        lo, hi = max(0, c - half), min(n, c + half)
        ecg[lo:hi] += beat_wave(t[lo:hi] - bt, kind, amp)
    ecg *= lead_gain
    # baseline wander: respiration (0.25 Hz) + slow drift (random walk, < 0.05 Hz)
    drift = np.cumsum(RNG.normal(0, 1, n))
    drift = sosfiltfilt(butter(1, 0.05, fs=FS, output='sos'), drift)
    drift = 0.25 * drift / (np.std(drift) + 1e-12)
    wander = 0.12 * np.sin(2 * np.pi * 0.25 * t + RNG.uniform(0, 6)) + drift
    mains = 0.015 * np.sin(2 * np.pi * 60.0 * t)
    return ecg + wander + mains + RNG.normal(0, noise_sd, n)


def corrupt(sig, t0, t1, rng_sd=0.6):
    """Electrode artefact on one lead: motion noise + two short flat dropouts."""
    s = sig.copy()
    i0, i1 = int(t0 * FS), int(t1 * FS)
    noise = sosfiltfilt(butter(2, [3, 30], 'bandpass', fs=FS, output='sos'), RNG.normal(0, rng_sd, i1 - i0))
    s[i0:i1] += noise * 2.0
    for start in RNG.uniform(t0 + 10, t1 - 10, 2):
        k = int(start * FS)
        s[k:k + int(1.5 * FS)] = s[k]
    return s


# ─── files ──────────────────────────────────────────────────────────────────
def write_edf(path, leads, t_from, t_to, rec_start):
    """Write samples [t_from, t_to) of the leads (mV) as an EDF file.
    Polarity is inverted, as on the Equivital (the pipeline sets INVERT_ECG)."""
    i0, i1 = int(round(t_from * FS)), int(round(t_to * FS))
    data = np.vstack([-leads[k][i0:i1] * 1e-3 for k in leads])  # mV -> V, inverted
    info = mne.create_info(list(leads), FS, ch_types='ecg')
    raw = mne.io.RawArray(data, info, verbose=False)
    raw.set_meas_date(rec_start + timedelta(seconds=t_from))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    mne.export.export_raw(path, raw, fmt='edf', overwrite=True, verbose=False)
    print(f'  wrote {os.path.relpath(path, OUT)}  ({t_to - t_from:.0f} s)')


def session_trials(state):
    """(name, start_s, duration_s) relative to the recording start."""
    s = state
    return [(f'{s}_5min_Sitting', 20, 300), (f'{s}_5min_Standing', 340, 300),
            (f'{s}_30s_T1', 660, 30), (f'{s}_30s_T2', 700, 30), (f'{s}_30s_T3', 740, 30)]


DURATION = 800.0  # s per session


def make_session(pid, state, rec_start, extra_events=(), bigem=None):
    blocks = [(0, 330, 0.85), (330, 650, 0.72), (650, DURATION, 0.80)]  # sitting / standing / pre-walk
    events = [(60, 'PVC'), (130, 'APC'), (200, 'MISS'), (265, 'PVC'),
              (400, 'PVC'), (470, 'NOISE'), (555, 'APC'),
              (701.0, 'PVC'), (712.5 if pid != 'P-001' else 705.0, 'MISS'), (755, 'PVC'), *extra_events]
    beats = build_beats(blocks, events)
    if bigem:
        beats = bigeminy(beats, *bigem)
    beats = beats[beats.time_s < DURATION - 0.5]
    leads = {'ECG 1': synth_ecg(beats, DURATION, 1.0), 'ECG 2': synth_ecg(beats, DURATION, 0.8, 0.025)}
    truth = beats.assign(clock=[rec_start.replace(tzinfo=None) + timedelta(seconds=float(x)) for x in beats.time_s])
    return leads, truth


def sync_rows(trials, rec_start):
    rows = []
    for name, t0, dur in trials:
        for ev, tt in (('start', t0), ('stop', t0 + dur)):
            ts = (rec_start + timedelta(seconds=tt)).strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]
            rows.append({'Timestamp': ts, 'Trial': name, 'Event': ev})
    return rows


def main():
    print(f'Generating synthetic recordings under {OUT}')
    all_sync = {}
    sheets = {}

    # ── P-001 ────────────────────────────────────────────────────────────────
    pid, edf = 'P-001', os.path.join(OUT, 'P-001', 'Equivital', 'EDF')
    off_start = SESSION_DATE.replace(hour=9)
    leads, truth = make_session(pid, 'OFF', off_start)
    write_edf(os.path.join(edf, 'P-001_OFF_State(1).edf'), leads, 0, 710, off_start)     # Equivital restart:
    write_edf(os.path.join(edf, 'P-001_OFF_State(2).edf'), leads, 720, DURATION, off_start)  # 10-s hole in T2 [D26]
    truth.to_csv(os.path.join(edf, 'P-001_OFF_truth.csv'), index=False)
    on_start = SESSION_DATE.replace(hour=11)
    leads, truth = make_session(pid, 'ON', on_start)
    write_edf(os.path.join(edf, 'PD_P001_ON_State_Rest(1).edf'), leads, 0, 326, on_start)  # sitting [D25]
    write_edf(os.path.join(edf, 'P-001_ON_State(1).edf'), leads, 330, DURATION, on_start)
    truth.to_csv(os.path.join(edf, 'P-001_ON_truth.csv'), index=False)
    all_sync[pid] = (sync_rows(session_trials('OFF'), off_start) + sync_rows(session_trials('ON'), on_start),
                     'Sync/SyncFile_Complete_Trial.csv')
    sheets[pid] = 'sub001_TrialRunSheet.csv'

    # ── P-015 ────────────────────────────────────────────────────────────────
    pid, edf = 'P-015', os.path.join(OUT, 'P-015', 'Equivital', 'EDF')
    off_start = SESSION_DATE.replace(hour=13)
    leads, truth = make_session(pid, 'OFF', off_start)
    leads['ECG 1'] = corrupt(leads['ECG 1'], 25, 315)  # poor ECG 1 during the sitting -> read ECG 2 [D27]
    write_edf(os.path.join(edf, 'P015_OFF(1).edf'), leads, 0, DURATION, off_start)
    truth.to_csv(os.path.join(edf, 'P015_OFF_truth.csv'), index=False)
    on_start = SESSION_DATE.replace(hour=15)
    leads, truth = make_session(pid, 'ON', on_start, bigem=(420, 480))  # 60-s bigeminy [D21/D28]
    write_edf(os.path.join(edf, 'P015_ON(1).edf'), leads, 0, DURATION, on_start)
    truth.to_csv(os.path.join(edf, 'P015_ON_truth.csv'), index=False)
    off_trials = session_trials('OFF') + [('ON_30s_T9', 775, 20)]  # ON trial inside the OFF files [D20]
    on_trials = session_trials('ON') + [('ON_IART_Baseline', 780, 2)]  # placeholder (< 5 s) [D25]
    all_sync[pid] = (sync_rows(off_trials, off_start) + sync_rows(on_trials, on_start),
                     'Sync/P-015_ses-01_sessionlog (2).csv')
    sheets[pid] = 'sub015_TrialRunSheet.csv'

    for pid, (rows, rel) in all_sync.items():
        path = os.path.join(OUT, pid, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        pd.DataFrame(rows).sort_values('Timestamp').to_csv(path, index=False)
        print(f'  wrote {os.path.relpath(path, OUT)}  ({len(rows) // 2} trials)')
        names = sorted({r['Trial'] for r in rows})
        notes = {'OFF_5min_Sitting': 'Use ECG 2 channel, NOT ECG 1 for the 5 minutes of quiet sitting'} \
            if pid == 'P-015' else {}
        pd.DataFrame({'Trial': names, 'Notes': [notes.get(n, '') for n in names]}).to_csv(
            os.path.join(OUT, pid, sheets[pid]), index=False)


if __name__ == '__main__':
    main()
