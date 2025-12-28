import os 
import sys
import numpy as np
import pandas as pd
from pandas.tseries.offsets import MonthEnd
from scipy.linalg import pinv, eigh
from scipy.stats import iqr
from scipy.sparse.linalg import eigsh
from scipy.stats import scoreatpercentile

# 2. Get current directory and split into parts
current_dir = os.getcwd()
parts = current_dir.split(os.sep)  # os.sep is the file separator ('/' or '\')

# 3. Reconstruct parent path (excluding last part)
parent_path = os.sep.join(parts[:-1])  # FORSE INUTILE

# 4. Add parent directory and its subfolders to Python path
sys.path.append(parent_path)

for root, dirs, files in os.walk(parent_path):
    sys.path.append(root)

sep = os.sep 

def main():
    
    sheet, IDf, trans, idt, method, q, algo = settings()

    if sheet == 'all':
       country = ['AT', 'BE', 'DE', 'EA', 'EL', 'ES', 'FR', 'IE', 'IT', 'NL', 'PT']
    else:
       country = [sheet]

     # Create output directory
    saveto = os.path.join(current_dir, f"data_TR{idt}")
    os.makedirs(saveto, exist_ok=True)

    for cc in range(len(country)):
        # ====================== #
        #   STEP 1: Load data    #
        # ====================== #
        sheet = country[cc]
        
        # Read Excel files using pandas
        try:
            # Load dataset (preserving column names exactly as in Excel)
            dataset = pd.read_excel(f"{sheet}data.xlsx", sheet_name='data')
            
            # Load info and convert to dictionary 
            info = pd.read_excel(f"{sheet}data.xlsx", sheet_name='info')
            spec = {col: info[col].values for col in info.columns}

            # Select transformations based on idt
            spec['TR'] = spec[f'TR{idt}']

            print(f"Successfully loaded data for {sheet}")
            
        except Exception as e:
            print(f"Error loading data for {sheet}: {str(e)}")
            continue  # Skip to next country if error occurs

        data = dataset.iloc[:, 1:].to_numpy(dtype=np.float64)  # Force float type
        titles = [str(x) for x in dataset.columns[1:]]  # Ensure strings
        dates = dataset.iloc[:, 0].to_numpy()  # Preserve original dtype
        
        # Optional date parsing if first column contains dates
        try:
            dates = pd.to_datetime(dates).to_numpy()  # Convert to datetime64
        except:
            pass  # Keep as-is if not date-like

        # ======================== #
        #   STEP 2: Set frequency  #
        # ======================== #
        if IDf == 'QM':
            # Case 1: quarterly aggregation (monthly to quarterly)
            data, dates, _ = aggregate(data, dates, spec) 
            dates = pd.to_datetime(dates).to_period('Q').to_timestamp().to_numpy()
        elif IDf == 'Q':
            # Case 2: only quarterly data
            locQ = np.array([freq == 'Q' for freq in spec['Frequency']])
            data = data[:, locQ]
            # Convert titles to numpy array for boolean indexing
            titles = np.array(titles)[locQ].tolist() if isinstance(titles, list) else titles[locQ]
        
           
            
            # Filter spec dictionary for quarterly series
            for key in ['TR', 'Class', 'Name', 'Aggregation', 'Frequency']:
                if key in spec:
                    spec[key] = np.asarray(spec[key])[locQ]
            
            data, dates, _ = aggregate(data, dates, spec)
            dates = pd.to_datetime(dates).to_period('Q').to_timestamp().to_numpy()

        elif IDf == 'M':
            # Case 3: only monthly data
            locM = np.array([freq == 'M' for freq in spec['Frequency']])
            data = data[:, locM]
            
            # Handle titles whether it's list or numpy array
            titles = np.array(titles)[locM].tolist() if isinstance(titles, list) else titles[locM]
            
            # Filter spec dictionary
            for key in ['TR', 'Class', 'Name', 'Aggregation', 'Frequency']:
                if key in spec:
                    spec[key] = np.asarray(spec[key])[locM]
        else:
            spec['agg'] = 0

        # ======================== #
        #  STEP 3: Transform Data  #
        # ======================== #
        spec['trans'] = trans
        xt = EA_transform(data, spec)

        # ===================== #
        #  STEP 4: Impute Data  #
        # ===================== #
        opts = {'maxiter': 1000, 'thresh': 1e-5, 'algo': algo}  # Set EM algorithm parameters

        # Convert dates to datetime if they're strings (assuming dates is a list/array)
        dates = pd.to_datetime(dates)  # Convert to datetime objects if needed

        # Find specific dates
        T19 = np.where(dates == pd.to_datetime('2019-10-01'))[0][0]  # Returns first index of Oct 1 2019
        T21 = np.where(dates == pd.to_datetime('2021-10-01'))[0][0]  # Returns first index of Oct 1 2021
        opts['T19'] = T19

        if method == -1:
            # Case -1: no imputation
            Xc = xt
            IDi = 'NA_'
            
        elif method == 0:
            # Case 0: impute only ragged edges
            opts['out'] = 0
            Xc = xt.copy()
            Xc, _ = EMimputation(xt, q, opts)
            IDi = ''
            
        elif method == 1:
            # Case 1: standard outlier imputation
            opts['out'] = 1
            Xc, pc = EMimputation(xt, q, opts)
            IDi = 'OA_'
            
        elif method == 2:
            # Case 2: impute real variables during Covid
            loc = [i for i, cls in enumerate(spec['Class']) if cls == 'R']
            Xnan = xt.copy()
            Xnan[T19+1:T21, loc] = np.nan
            Xc, pc = EMimputation(Xnan, q, opts)
            IDi = 'COV_'
            
        # Create and save timetable
        TT = pd.DataFrame(Xc, index=dates, columns=spec['Name'])
        TT.index.name = 'Date'
        
        # Define file name
        sIDi = f"_{IDi}" if IDi else IDi
        nsave = os.path.join(saveto, f"{sheet}data{IDf}{sIDi}_TR{idt}.xlsx")        
        
        # Save to Excel
        TT.to_excel(nsave, sheet_name=sheet)
        
        print(f"Processed and saved {sheet}")

    return Xc



def settings():
    """
    Command-line interface for model settings.
    
    Returns:
        sheet (str): Country identifier ('all' or specific)
        IDf (str): Frequency ('M', 'Q', or 'QM')
        trans (str): Transformation method ('light' or 'heavy')
        method (int/str): Imputation method (number or empty string)
        q (int): Number of factors for imputation
    """
    # Ask for default options
    use_default = input("Use default options? (y/n): ").strip().lower()
    
    if use_default == 'y':
        sheet = 'all'
        IDf = 'QM'
        trans = 'light'
        method = 0
        q = 99
        algo = 'SW'
    else:
        # Country selection
        sheet = input("Select country [default: all]: ").strip()
        if not sheet:
            print("No country specified, downloading data for all countries")
            sheet = 'all'
        
        # Frequency selection
        IDf = input("Select frequency (QM/Q/M/MF) [default: QM]: ").strip().upper()
        if IDf not in ('QM', 'Q', 'M', 'MF'):
            print("Invalid frequency, using <QM>")
            IDf = 'QM'
        
        # Transformation selection
        trans = input("Select transformation (light/heavy/BLT) [default: light]: ").strip().lower()
        if trans not in ('light', 'heavy', 'blt'):
            print("Invalid transformation, using <light>")
            trans = 'light'
        
        # Imputation selection
        imp = input("Impute missing values/outliers/Covid period? (y/n) [default: y]: ").strip().lower()
        if imp != 'n':
            method = input("Imputation method [default: 0]: ").strip()
            method = int(method) if method and method.isdigit() else 0
            
            q = input("Number of factors [default: 99]: ").strip()
            q = int(q) if q and q.isdigit() else 99
            
            # Algorithm selection
            algo = input("Select imputation algorithm (SW/BM) [default: SW]: ").strip().upper()
            if algo not in ('SW', 'BM'):
                print("Invalid algorithm, using 'SW'")
                algo = 'SW'
        else:
            method = -1
            q = 99
            algo = ''
    
    # Set transformation identifier
    if trans == 'light':
        idt = 2
    elif trans == 'heavy':
        idt = 1
    elif trans == 'blt':
        idt = 3
    
    # Mixed-frequency imputation warning
    if IDf == 'MF' and method > -1:
        print("Warning: Imputation of missing values not available for mixed-frequency data, setting <method> to -1")
        method = -1
    
    # Print summary
    print(f"Country: {sheet}; Frequency: {IDf}; trans: {trans}; method: {method}; q: {q}; algorithm: {algo}.")
    
    return sheet, IDf, trans, idt, method, q, algo

################################################################################

def aggregate(data, dates, spec):
    """
    Aggregates mixed-frequency (monthly-quarterly) data to quarterly frequency
    
    Parameters:
    -----------
    data : ndarray
        Input data (TxN) with NaNs
    dates : array-like
        Datetime objects (Tx1)
    spec : dict
        Contains:
        - Frequency: List of 'Q'/'M' for each variable
        - Aggregation: List of 1 (mean) or 2 (sum) for each variable
        - Name: Optional list of variable names
    
    Returns:
    --------
    dataQ : ndarray
        Quarterly aggregated data
    datesQ : ndarray
        Quarterly dates
    datasetQ : DataFrame
        DataFrame with dates as index
    """
    
    # Input validation
    n_vars = data.shape[1]

    if len(spec['Frequency']) != n_vars or len(spec['Aggregation']) != n_vars:
        raise ValueError('Identifiers must match data dimension')
    
    # Convert to DataFrame
    dates = pd.to_datetime(dates.flatten())
    temp = data.copy()
    tm = dates[len(dates)//2]
    
    # Process monthly variables
    for i in [i for i, freq in enumerate(spec['Frequency']) if freq == 'M']:
        idxNaN = np.isnan(data[:, i])
        valid_dates = dates[~idxNaN]
        
        if len(valid_dates) == 0:
            continue
            
        t0, t1 = min(valid_dates), max(valid_dates)
        
        # Handle single observation case
        if t0 == t1:
            if t0 < tm:
                t1 = None
            else:
                t0 = None
        
        # Check start of quarter
        if t0 is not None:
            tp = t0 + MonthEnd(1)
            tp1 = tp + MonthEnd(1)
            if (tp.quarter != t0.quarter) or (tp1.quarter != t0.quarter):
                temp[dates == t0, i] = np.nan
        
        # Check end of quarter
        if t1 is not None:
            tp = t1 - MonthEnd(1)
            tp1 = tp - MonthEnd(1)
            if (tp.quarter != t1.quarter) or (tp1.quarter != t1.quarter):
                temp[dates == t1, i] = np.nan
    
    # Create DataFrame
    var_names = spec.get('Name', [str(i) for i in range(n_vars)])
    df = pd.DataFrame(temp, index=dates, columns=var_names)
    
    # Identify flow variables (sum aggregation)
    if 'M' in spec['Frequency']:
        flow_vars = [var_names[i] for i, (freq, agg) in enumerate(zip(spec['Frequency'], spec['Aggregation'])) 
                     if freq == 'M' and agg == 2]
    else:
        flow_vars = []
    
    # Resample to quarterly
    df_q = df.resample('QS').mean()  # Default is mean (stock variables)
    # Handle flow variables (sum aggregation)
    for var in flow_vars:
        df_q[var] = df[var].resample('QS').sum()

    # 2. Shift labels to middle month (+2 months)
    df_q.index = df_q.index + pd.DateOffset(months=2)  # Move to 2nd month of quarter
    
    # Prepare outputs
    dataQ = df_q.to_numpy()
    datesQ = df_q.index.to_numpy()
    datasetQ = df_q.reset_index()
    
    return dataQ, datesQ, datasetQ

###############################################################################

def mdiff(X, s=1, k=1):
    """
    Performs k-th differences with a mixed-frequency panel of data
    
    Parameters:
    -----------
    X : ndarray
        Matrix of data with NaNs (TxN)
    s : int
        Period length (2 for quarters, 11 for years)
    k : int
        Order of differencing (default: 1)
    
    Returns:
    --------
    Y : ndarray
        Matrix of differenced data (TxN)
    """
    
    Y = np.full_like(X, np.nan)  # Initialize output with NaNs
    
    # Handle 1st order differences
    dX = X[s:,:] - X[:-s,:]  # 1st differences between periods s
    
    # Handle higher-order differences
    for _ in range(k-1):
        dX = dX[s:,:] - dX[:-s,:]
    
    # Fill output matrix
    Y[s*k:,:] = dX
    
    return Y

import numpy as np

###############################################################################

def remove_outliers(X, c=10):
    """
    Identifies and replaces outliers with NaNs based on the median and interquartile range.
    
    Parameters:
    -----------
    X : ndarray
        Input data matrix (TxN)
    c : float
        Threshold multiplier for outlier detection (default=10)
    
    Returns:
    --------
    X_cleaned : ndarray
        Data with outliers replaced by NaNs
    out : ndarray (bool)
        Boolean mask of outlier positions
    n : ndarray
        Number of outliers per column
    """
    # Calculate median and IQR for each column (ignoring NaNs)
    mX = np.nanmedian(X, axis=0)
    q75, q25 = np.nanpercentile(X, [75, 25], axis=0, method='linear')
    iqrX = q75 - q25 # most similar to matlab

    # Identify outliers (abs(X - median) > c * IQR)
    out = np.abs(X - mX) > c * iqrX
    
    # Replace outliers with NaN
    X_clean = X.copy()
    X_clean[out] = np.nan

    # Control for series with many zeros (>20% outliers)
    out_ratio = np.sum(out, axis=0) / X.shape[0]
    loc = out_ratio > 0.2
    out[:, loc] = False
    
    # Count outliers per column
    n = np.sum(out, axis=0)
    
    return X_clean, out, n

###############################################################################

def EA_transform(X, spec, c=1):
    """
    Applies transformations to variables based on spec['TR'] codes
    
    Parameters:
    -----------
    X : ndarray (TxN)
        Input data matrix
    spec : dict
        - TR: Transformation codes (1-6) for each variable
        - trans: 'light' or 'heavy' (for fractional TR handling)
        - Frequency: List of 'Q'/'M' for each variable
        - agg: 1 (balanced) or 0 (mixed frequencies)
    c : float
        Scaling factor (default=1)
    
    Returns:
    --------
    Xt : ndarray (TxN)
        Transformed data
    """
    # Initialize
    spec.setdefault('agg', 1)
    X = np.array(X, dtype=float)
    Xt = np.full_like(X, np.nan)
  
    # Validate
    if len(spec['TR']) != X.shape[1]:
        if len(spec['TR']) == X.shape[0]: 
            X = X.T
        else:
            raise ValueError("TR length must match X columns")
    
    # Check for negative values where logs are needed
    log_trs = [1, 2, 3]
    neg_mask = (X < 0) & np.isin(spec['TR'], log_trs)
    
    if np.any(neg_mask):
        bad_cols = np.where(np.any(neg_mask, axis=0))[0]
        print(f"Warning: Variables in positions {bad_cols} contain negative values")
        print("Logs are not admitted. Transforming without logs")
        
        # Replace log transformations with non-log equivalents
        for col in bad_cols:
            if spec['TR'][col] == 1:
                spec['TR'][col] = 4
            elif spec['TR'][col] == 2:
                spec['TR'][col] = 5
            elif spec['TR'][col] == 3:
                spec['TR'][col] = 6
    
    TR = spec['TR']
    # Apply transformations
   # Apply transformations
    if spec['agg'] == 1:  # Balanced panel
        masks = {
            1: TR == 1,  # c*log(x)
            2: TR == 2,  # Δlog(c*x)
            3: TR == 3,  # ΔΔlog(c*x)
            4: TR == 4,  # x (no transform)
            5: TR == 5,  # Δx
            6: TR == 6   # ΔΔx
        }
        
        Xt[:, masks[1]] = c * np.log(X[:, masks[1]])
        Xt[1:, masks[2]] = np.diff(c * np.log(X[:, masks[2]]), axis=0)
        Xt[2:, masks[3]] = np.diff(c * np.log(X[:, masks[3]]), n=2, axis=0)
        Xt[:, masks[4]] = X[:, masks[4]]
        Xt[1:, masks[5]] = np.diff(X[:, masks[5]], axis=0)
        Xt[2:, masks[6]] = np.diff(X[:, masks[6]], n=2, axis=0)
        
    else:  # Mixed frequencies
        is_m = np.array(spec['Frequency']) == 'M'
        is_q = ~is_m
        
        # TR=1 (logs)
        mask = TR == 1
        Xt[:, mask] = c * np.log(X[:, mask])
        
        # TR=2 (Δlog)
        Xt[1:, is_m & (TR == 2)] = np.diff(c * np.log(X[:, is_m & (TR == 2)]), axis=0)
        Xt[:, is_q & (TR == 2)] = mdiff(c * np.log(X[:, is_q & (TR == 2)]), 3, 1)
        
        # TR=3 (ΔΔlog)
        Xt[2:, is_m & (TR == 3)] = np.diff(c * np.log(X[:, is_m & (TR == 3)]), n=2, axis=0)
        Xt[:, is_q & (TR == 3)] = mdiff(c * np.log(X[:, is_q & (TR == 3)]), 3, 2)
        
        # TR=4 (no transform)
        Xt[:, TR == 4] = X[:, TR == 4]
        
        # TR=5 (Δx)
        Xt[1:, is_m & (TR == 5)] = np.diff(X[:, is_m & (TR == 5)], axis=0)
        Xt[:, is_q & (TR == 5)] = mdiff(X[:, is_q & (TR == 5)], 3, 1)
        
        # TR=6 (ΔΔx)
        Xt[2:, is_m & (TR == 6)] = np.diff(X[:, is_m & (TR == 6)], n=2, axis=0)
        Xt[:, is_q & (TR == 6)] = mdiff(X[:, is_q & (TR == 6)], 3, 2)
    
    return Xt

###############################################################################

def EMimputation(X, q0, opts=None):
    """
    Python implementation of EM imputation using factor models
    
    Parameters:
    -----------
    X : ndarray
        Input data matrix (TxN) with missing values as NaN
    q0 : int
        Number of factors (99 for Bai & Ng IC selection)
    opts : dict, optional
        Options dictionary with:
        - maxiter: Maximum iterations (default 500)
        - thresh: Convergence threshold (default 1e-6)
        - out: Outlier removal flag (default 0)
    
    Returns:
    --------
    X : ndarray
        Data with imputed values
    pc : dict
        Factor model results containing:
        - F: Factors (Txq)
        - C: Loadings (Nxq)
        - chi: Common component (TxN)
        - R: Idiosyncratic variances (N,)
        - d: Eigenvalues (q,)
    """
    # Set default options
    if opts is None:
        opts = {}
    opts.setdefault('maxiter', 1000)
    opts.setdefault('thresh', 1e-5)
    opts.setdefault('out', 0)

    if opts is None:   
        print(f"Maximum number of iterations set to {opts['maxiter']}")
        print(f"Threshold for EM objective set to {opts['thresh']:10f}")
    
    # Outlier removal
    if opts['out'] != 0:
        X, _, _ = remove_outliers(X)
    
    # Initial imputation with means 
    indNaN = np.isnan(X)
    mx = np.tile(np.nanmean(X, axis=0), (X.shape[0], 1))  # data mean, no NaNs
    sx = np.tile(np.nanstd(X, axis=0), (X.shape[0], 1))   # data std, no NaNs
    X[indNaN] = mx[indNaN]
    
    # Standardize
    mx = np.tile(np.mean(X, axis=0), (X.shape[0], 1))
    sx = np.tile(np.std(X, axis=0, ddof=1), (X.shape[0], 1))
    Xs = (X - mx) / sx
    
    # Determine number of factors
    q = BaiNg(Xs, 15, 2) if q0 == 99 else q0
    
    # Initial PCA
    pc = princfact(Xs, q)
    chi0 = pc['chi']
    
    EM = {}
    EM['chi'] = pc['chi']
    
    # EM algorithm
    err = 999
    j = 0
    
    while err > opts['thresh'] and j < opts['maxiter']:
        if j % 10 == 0:
            print(f'Running iteration {j}: error is {err:10f}')
        
        # Imputation step
        X[indNaN] = EM['chi'][indNaN] * sx[indNaN] + mx[indNaN]
        
        # Re-standardize
        mx = np.tile(np.mean(X, axis=0), (X.shape[0], 1))
        sx = np.tile(np.std(X, axis=0, ddof=1), (X.shape[0], 1))
        Xs = (X - mx) / sx
        
        # Update number of factors if needed
        q = BaiNg(Xs, 20, 2) if q0 == 99 else q0
        
        # PCA step
        EM = princfact(Xs, q)
        
        # Compute convergence criterion
        dchi = EM['chi'] - chi0
        err = np.sum(dchi**2) / np.sum(chi0**2)
        
        # Update for next iteration
        chi0 = EM['chi']
        j += 1
        
        if err < opts['thresh'] and j < opts['maxiter']:
            print(f'EM converged after {j} iterations')
    
    return X, EM

###############################################################################
def princfact(X, q, method=2, st=1):
    """
    Estimates factors via principal components with NaN handling
    
    Parameters:
    -----------
    X : ndarray
        Input data matrix (TxN)
    q : int
        Number of factors to extract
    method : int
        Normalization method (0-3, default=2)
    st : int
        Standardization flag (0/1, default=0)
    
    Returns:
    --------
    pc : dict
        Contains:
        - F: Factors (Txq)
        - C: Loadings (Nxq)
        - chi: Common component (TxN)
        - R: Idiosyncratic variances (N,)
        - d: Eigenvalues (q,)
    """
    T, N = X.shape
    
    # Handle NaNs
    if np.isnan(X).any():
        X = np.where(np.isnan(X), np.nanmean(X), X)
        print('Detected NaNs in original data matrix: replaced with unconditional mean of non-NaNs values')
    
    # Standardize if requested
    if st == 1:
        X = (X - np.mean(X, axis=0)) / np.std(X, axis=0, ddof=1)
    
    # Validate method
    if method not in [0, 1, 2, 3]:
        print('The value for <method> is outside the bounds, setting <method> to 2')
        method = 2
    
    # Compute covariance matrix
    cov_mat = np.cov(X.T, bias = True)
    
    # Factor extraction
    if method != 3:
        # Standard NxN covariance approach
        d, v = eigsh(cov_mat, k = q, which = 'LM')  # Get q largest eigenvalues/vectors
        d = np.flip(d)
        v = np.flip(v, axis=1)
        
        # Sign flip for consistency
        if np.sum(v) < 0:
            v = -v
            
        if method == 0:
            f = X @ v
            C = v
        elif method == 1:
            f = X @ v @ np.diag(1/np.sqrt(d))
            C = v @ np.diag(np.sqrt(d))
        elif method == 2:
            f = X @ v / np.sqrt(N)
            C = v * np.sqrt(N)
    else:
        # TxT covariance approach (method 3)
        cov_mat_t = np.cov(X, bias = True)
        d, v = eigsh(cov_mat_t, k = q, which='LM')  # Get q largest eigenvalues/vectors
        d = np.flip(d)
        v = np.flip(v, axis=1)
        
        if np.sum(v) < 0:
            v = -v
            
        f = v * np.sqrt(T)
        C = (f.T @ X).T / T
    
    # Common component and residuals
    chi = f @ C.T
    R = np.diag(np.diag(np.cov(X - chi, rowvar=False)))
    
    return {
        'F': f,
        'C': C,
        'chi': chi,
        'R': R,
        'd': d
    }

###############################################################################
def BaiNg(X, qmax=None, method=2):
    """
    Selects number of factors using Bai & Ng (2002) information criteria
    
    Parameters:
    -----------
    X : ndarray
        Input data matrix (TxN)
    qmax : int
        Maximum number of factors to test (default=sqrt(N))
    method : int
        Criterion method (1-3, default=2)
    
    Returns:
    --------
    q : int
        Number of factors selected
    """
    T, N = X.shape
    
    # Set defaults
    if qmax is None:
        qmax = int(np.sqrt(N))
    
    # Handle NaNs
    if np.isnan(X).any():
        X = np.where(np.isnan(X), np.nanmean(X), X)
        print('Detected NaNs in original data matrix: replaced with unconditional mean of non-NaNs values')
    
    # Standardize data
    X = (X - np.ones((T, 1)) @ (np.sum(X, axis=0).reshape(1, -1) / T))/ np.std(X, axis=0, ddof=1)
    
    # Initialize criterion
    CT = np.zeros(qmax)
    
    # Select penalty
    if method == 1:
        CT = np.log(N*T/(N+T)) * np.arange(1, qmax+1) * (N+T)/(N*T)
    elif method == 2:
        CT = ((N+T)/(N*T)) * np.log(min(N, T)) * np.arange(1, qmax+1)
    elif method == 3:
        CT = np.arange(1, qmax+1) * np.log(min(N, T))/min(N, T)
    else:
        raise ValueError("Method must be 1, 2, or 3")
    
    # Principal components
    if T < N:
        # Use TxT covariance matrix
        cov_mat = np.cov(X)
        d, v = eigsh(cov_mat, k=qmax, which='LM')  # Get qmax largest eigenvalues/vectors
        d = np.flip(d)
        v = np.flip(v, axis=1)
        Fhat = np.sqrt(T) * v
        Lhat = X.T @ Fhat / T
    else:
        # Use NxN covariance matrix
        cov_mat = np.cov(X.T, bias = True)
        d, v = eigsh(cov_mat,k=qmax, which='LM')  # Get qmax largest eigenvalues/vectors
        d = np.flip(d)
        v = np.flip(v, axis=1)
        Lhat = np.sqrt(N) * v
        Fhat = X @ Lhat / N

    # Information criterion
    Sigma = np.zeros(qmax+1)
    IC1 = np.zeros(qmax+1)
    
    for qq in range(qmax, 0, -1):
        chat = Fhat[:, :qq] @ Lhat[:, :qq].T
        ehat = X - chat
        Sigma[qq-1] = np.mean(np.sum(ehat * ehat, axis=0) / T)
        IC1[qq-1] = np.log(Sigma[qq-1]) + CT[qq-1]
    
    # Case with no factors
    Sigma[qmax] = np.mean(np.sum(X * X, axis=0) / T)
    IC1[qmax] = np.log(Sigma[qmax])

    # Select number of factors
    q = np.argmin(IC1)+1
    if q > qmax:
        q = 0

    return q

###############################################################################
def kalman(Y, pars):
    """
    Python implementation of MATLAB's Kalman Filter and Smoother
    Maintains identical input/output structure and variable names
    """
    # Unpack parameters
    A = pars['A']
    C = pars['C']
    R = pars['R']
    Q = pars['Q']
    
    # Initialize parameters
    Z0 = pars.get('Z00', np.zeros((C.shape[1], 1)))
    P0 = pars.get('P00', 10 * np.eye(C.shape[1]))
    mu = pars.get('mu', np.zeros((C.shape[1], 1)))
    beta = pars.get('beta', np.zeros((1, 2)))
    
    ns = C.shape[1]  # Number of states
    Y = Y.T  
    N, T = Y.shape
    
    # Initialize output structure 
    KF = {
        'Zttm': np.nan * np.zeros((ns, T)),
        'Pttm': np.nan * np.zeros((ns, ns, T)),
        'Ztt': np.nan * np.zeros((ns, T+1)),
        'Ptt': np.nan * np.zeros((ns, ns, T+1)),
        'Kt': np.zeros((ns, N, T)),
        'vt': np.zeros((N, T)),
        'loglik': 0,
        'Z00': Z0,
        'P00': P0
    }
    
    KF['Ztt'][:, [0]] = Z0
    KF['Ptt'][:, :, 0] = P0
    mu = np.hstack([np.zeros((ns, 1)), np.tile(mu, (1, T-1))])
    
    # Kalman Filter %%
    for t in range(T):
        # Prediction Step
        KF['Zttm'][:, t] = mu[:, t] + A @ KF['Ztt'][:, t]
        KF['Pttm'][:, :, t] = A @ KF['Ptt'][:, :, t] @ A.T + Q
        KF['Pttm'][:, :, t] = 0.5 * (KF['Pttm'][:, :, t] + KF['Pttm'][:, :, t].T)
        
        # Update Step
        idx = ~np.isnan(Y[:, t])
        yt = Y[idx, t]
        Ct = C[idx, :]
        Rt = R[np.ix_(idx, idx)]
        bt = beta[0, idx] if np.any(beta) else beta
        
        if len(yt) == 0:
            KF['Ztt'][:, t+1] = KF['Zttm'][:, t]
            KF['Ptt'][:, :, t+1] = KF['Pttm'][:, :, t]
        else:
            Kt_part = KF['Pttm'][:, :, t] @ Ct.T @ np.linalg.pinv(Ct @ KF['Pttm'][:, :, t] @ Ct.T + Rt)
            KF['Kt'][:, idx, t] = Kt_part
            KF['vt'][idx, t] = yt - Ct @ KF['Zttm'][:, [t]] - bt @ np.array([[1], [t+1]])
            
            KF['Ztt'][:, t+1] = KF['Zttm'][:, t] + KF['Kt'][:, :, t] @ KF['vt'][:, [t]]
            KF['Ptt'][:, :, t+1] = KF['Pttm'][:, :, t] - Kt_part @ Ct @ KF['Pttm'][:, :, t]
            KF['Ptt'][:, :, t+1] = 0.5 * (KF['Ptt'][:, :, t+1] + KF['Ptt'][:, :, t+1].T)
            
            # Log-likelihood calculation
            innov_cov = Ct @ KF['Pttm'][:, :, t] @ Ct.T + Rt
            KF['loglik'] += 0.5 * (np.log(np.linalg.det(np.linalg.pinv(innov_cov))) - 
                                  KF['vt'][idx, t].T @ np.linalg.pinv(innov_cov) @ KF['vt'][idx, t])
    
    #  Kalman Smoother %%
    KF['ZtT'] = np.zeros((ns, T))
    KF['PtT'] = np.zeros((ns, ns, T))
    KF['PtTm'] = np.zeros((ns, ns, T))
    
    KF['ZtT'][:, T-1] = KF['Ztt'][:, T]
    KF['PtT'][:, :, T-1] = KF['Ptt'][:, :, T]
    KF['PtTm'][:, :, T-1] = (np.eye(ns) - KF['Kt'][:, :, T-1] @ C) @ A @ KF['Ptt'][:, :, T-1]
    
    Pts = np.zeros((ns, ns, T))
    Pts[:, :, T-1] = KF['Ptt'][:, :, T-1] @ A.T @ pinv(KF['Pttm'][:, :, T-1])
    
    for t in range(T-1, 0, -1):
        KF['ZtT'][:, t-1] = (KF['Ztt'][:, t] + 
                             Pts[:, :, t] @ (KF['ZtT'][:, t] - KF['Zttm'][:, t]))
        
        KF['PtT'][:, :, t-1] = (KF['Ptt'][:, :, t] + 
                                Pts[:, :, t] @ (KF['PtT'][:, :, t] - KF['Pttm'][:, :, t]) @ Pts[:, :, t].T)
        
        if t > 1:
            Pts[:, :, t-1] = KF['Ptt'][:, :, t-1] @ A.T @ pinv(KF['Pttm'][:, :, t-1])
            KF['PtTm'][:, :, t-1] = (KF['Ptt'][:, :, t] @ Pts[:, :, t-1].T + 
                                     Pts[:, :, t] @ (KF['PtTm'][:, :, t] - A @ KF['Ptt'][:, :, t]) @ Pts[:, :, t-1].T)
    
    KF['ZtT'] = KF['ZtT'].T  
    
    return KF

###################################################################################


if __name__ == '__main__':
    Xc = main()












