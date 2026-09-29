import os
import re
from glob import glob

import mne
import numpy as np
import pandas as pd
import fooof
from neurodsp.spectral import compute_spectrum

# CONFIGURACIÓN
# Todas las rutas son relativas a la raíz del repositorio, por lo que el
# script funciona desde cualquier directorio de trabajo.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

main_dir = os.path.join(REPO_ROOT, "Data", "Prepro_data")
country = "Chile"
subject_class = "CN_02"

# Base conductual (record_id, PSS_DIRt, GHQ_DIRt, MAIA_Autorregulacion_DIRd)
behavioral_data_path = os.path.join(REPO_ROOT, "Data", "Behavioral_data.csv")

# Resultados intermedios (exponente por sujeto y lista de errores)
output_dir = os.path.join(
    REPO_ROOT, "Results_Aperiodic_2026", country, subject_class
)
os.makedirs(output_dir, exist_ok=True)

global_exponent_path = os.path.join(
    output_dir, f"{subject_class}_Exponent_Global.csv"
)

# Salida final: base de análisis que usan los scripts de R
merged_output_path = os.path.join(REPO_ROOT, "Analisis", "df.csv")

#
def normalize_record_id(value):
    """Convierte identificadores como subject_id_001 a Subject_ID001."""
    value = str(value).strip()
    match = re.search(r"(\d+)$", value)
    if match:
        return f"Subject_ID{match.group(1).zfill(3)}"
    return value


subjects_dir = os.path.join(main_dir, country, subject_class)
subjects_list = sorted([
    folder
    for folder in os.listdir(subjects_dir)
    if os.path.isdir(os.path.join(subjects_dir, folder))
])

print(f"Subjects found: {len(subjects_list)}")

results = []
error_list = []

for subj in subjects_list:
    record_id = normalize_record_id(subj)
    print(f"\n→ Processing subject: {record_id}")

    search_pattern = os.path.join(
        subjects_dir, subj, "eeg", "*.set"
    )
    files = sorted(glob(search_pattern))

    if not files:
        print("Error: No .set file found.")
        error_list.append(record_id)
        continue

    try:
        eeg_raw = mne.io.read_raw_eeglab(
            files[0], preload=True, verbose=False
        )
        fs = eeg_raw.info["sfreq"]
        eeg_data = eeg_raw.get_data()

        # PSD mediante Welch: ventanas de 1 segundo y 50% de solapamiento.
        freqs, psd = compute_spectrum(
            eeg_data,
            fs,
            method="welch",
            nperseg=int(fs),
            noverlap=int(fs / 2),
        )

        # Mismos parámetros FOOOF del script original.
        fg = fooof.FOOOFGroup(
            peak_width_limits=[1, 8],
            max_n_peaks=6,
            min_peak_height=0.2,
            peak_threshold=2.0,
            aperiodic_mode="fixed",
        )
        fg.fit(
            freqs=freqs,
            power_spectra=psd,
            freq_range=[2, 40],
        )

        # En modo fixed: columna 0 = offset; columna 1 = exponent.
        aperiodic_params = fg.get_params("aperiodic_params")
        channel_exponents = aperiodic_params[:, 1]

        # Exponente global: promedio de los exponentes de todos los canales.
        exponent_global = np.mean(channel_exponents)

        results.append({
            "record_id": record_id,
            "Exponent__Global": exponent_global,
        })

    except Exception as error:
        print(f"Error processing EEG: {error}")
        error_list.append(record_id)

global_exponent = pd.DataFrame(
    results, columns=["record_id", "Exponent__Global"]
)
global_exponent.to_csv(
    global_exponent_path, index=False, encoding="utf-8-sig"
)

if error_list:
    error_path = os.path.join(
        output_dir, f"{subject_class}_Errors.txt"
    )
    with open(error_path, "w", encoding="utf-8") as file:
        file.write("\n".join(error_list))
    print(f"Errors saved to: {error_path}")

print(f"Global exponent saved to: {global_exponent_path}")

#MERGE CON LA BASE CONDUCTUAL
data = pd.read_csv(behavioral_data_path)

if "record_id" not in data.columns:
    raise KeyError("La base conductual no contiene la columna 'record_id'.")

data["record_id"] = data["record_id"].map(normalize_record_id)

if data["record_id"].duplicated().any():
    duplicated = data.loc[
        data["record_id"].duplicated(keep=False), "record_id"
    ].unique()
    raise ValueError(
        f"Existen record_id duplicados en la base conductual: {duplicated.tolist()}"
    )

if global_exponent["record_id"].duplicated().any():
    raise ValueError("Existen record_id duplicados en los resultados EEG.")

if "Exponent__Global" in data.columns:
    raise ValueError(
        "La base conductual ya contiene una columna llamada 'Exponent__Global'."
    )

merged_data = data.merge(
    global_exponent,
    on="record_id",
    how="left",
    validate="one_to_one",
)

# Mismo orden de columnas que usan los scripts de análisis
first_cols = ["record_id", "Exponent__Global"]
merged_data = merged_data[
    first_cols + [c for c in merged_data.columns if c not in first_cols]
]

merged_data.to_csv(merged_output_path, index=False)

n_matched = merged_data["Exponent__Global"].notna().sum()
print(f"Rows in behavioral data: {len(data)}")
print(f"Rows with global exponent: {n_matched}")
print(f"Rows without global exponent: {len(data) - n_matched}")
print(f"Merged data saved to: {merged_output_path}")