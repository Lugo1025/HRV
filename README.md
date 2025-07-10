
ICU Survival & RR Interval Analyzer
===================================

This application provides an integrated GUI for clinicians and researchers to perform ICU survival prediction using clinical variables and heart rate variability (HRV) analysis from RR interval data. It combines machine learning models (ANN and Logistic Regression) and standard HRV analysis (Time, Frequency, and Nonlinear domains) into a user-friendly interface built with Tkinter.

Features
--------

- ICU Survival Prediction:
  - Predicts patient survival (Alive/Dead) using pre-trained ANN models based on selected HRV features.
  - Supports four model types:
    - Correlation (Top 15)
    - Correlation (Top 5)
    - Random Forest (Top 15)
    - Mutual Information (Top 15)
  - Provides logistic regression probability for interpretability.

- RR Interval Analysis:
  - Reads RR interval data from `.txt` files (space- or tab-delimited).
  - Computes time-domain, frequency-domain, and nonlinear HRV metrics using the `pyhrv` library.
  - Saves analysis results to a new `.txt` file per patient.
  - Visualizes RR time series in the GUI with interactive plots.

File Structure
--------------

models/
├── <model_name>_model.h5         # Trained ANN model
├── <model_name>_scaler.pkl       # Scaler used for feature normalization
└── <model_name>_features.txt     # Feature names used in model

gui general.py                    # Main Python GUI application

Requirements
------------

Make sure you have the following Python packages installed:

    pip install numpy pandas matplotlib tensorflow joblib pyhrv

If you're running on a Linux system without a display, set up your environment with:

    sudo apt-get install python3-tk

How to Run
----------

    python "gui general.py"

This will launch a GUI with two tabs:

1. ICU Prediction: 
   - Upload a `.txt` file containing clinical features (one line, whitespace-separated).
   - Select the model type.
   - Enter patient ID and name.
   - Click Analyze to predict survival and view ANN and logistic regression results.

2. RR Interval Analyzer:
   - Upload a `.txt` file with RR intervals and timestamps.
   - The app will extract HRV features, save them to a file, and plot the RR time series.

Input File Formats
------------------

### Clinical File (.txt)

A plain text file with a single row of numerical features. The number and order of features must match the selected model's expected input (see corresponding `features.txt`).

Example:

    793.8 75.2 56.4 92.1 123.4 ...

### RR Interval File (.txt)

A tab- or space-separated file with at least two columns: one for time and one for RR intervals.

Example:

    Time(ms)    RR(ms)
    0.00        812
    0.81        793
    1.60        785

Output
------

- RR analysis result is saved as:

      <PatientID>_<PatientName>_rr.txt

- ANN Prediction and Logistic Regression Probability are displayed in the GUI after analysis.

Customization
-------------

You can extend or retrain models by placing new `.h5`, `.pkl`, and `.txt` files in the `models/` directory and updating the `model_options` and `logistic_configs` dictionaries in the code.
