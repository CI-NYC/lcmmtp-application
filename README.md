# lcmmtp-application

Code to analyze the direct and indirect effects through acute kidney injury (AKI) of a delay in intubation policy on COVID-19 mortality using the lcmmtp R package.

## Data structure creation

1. First, we filtered for cohort eligibility criteria (e.g. age >= 18, no chronic kidney disease) in `scripts/0_create_cohort.R`.

2. Then, we selected for the relevant date-times for creating the wide, discrete time window data structure in later steps in `scripts/1-create_date_times.R`. These included:

- `t1_start` : start time of study, e.g. beginning of $L_1$ variables.

- `A_time` : time of exposure, intubation, if applicable

- `M_time` : time of mediator, AKI, if applicable

- `Y_time` : time out outcome, death, if applicable

- `Cens_time` : time out censoring event, discharge or transfer, if applicable

At this step we also filtered out <0.6% of patients with non-sensical date-times due to errors in data collection.

3. Next, we created variables `window`, `l_start`, `l_end`, `z_start`, and `z_end` to indicate what the date-time should be for each time window $t$ and random variable $L$ and $Z$. To enforce temporality in the time-discretized data structure, intervals depended on whether the patient was intubated (before AKI, if applicable), met AKI criteria (before intubation, if applicable), or was neither intubated nor met AKI criteria during the study follow up.

### Time intervals for patients who **never** were intubated nor met criteria for AKI:
