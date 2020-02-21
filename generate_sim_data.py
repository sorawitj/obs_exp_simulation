import pandas as pd
from utils import *
from statsmodels.distributions.empirical_distribution import ECDF


# Adapt from the R code in https://sites.google.com/view/ACIC2019DataChallenge/data-challenge
# Simulating observational and clinal trial data based on credit card fraud dataset
# The dataset consists of 30000 observations with 24 features
# binary outcome and treatment ('default payment next month' -> outcome, 'marriage' -> treatment)
# treatment and outcome are a complex function of measured covariates
# treatment effect is heterogeneous
# covariate distribution in clinical trial data is different from the one in observational data


def gen_population():
    # Yeh, I. C., & Lien, C. H. (2009). The comparisons of data mining techniques for the predictive accuracy of probability of default of credit card clients. Expert Systems with Applications, 36(2), 2473-2480.
    # https://archive.ics.uci.edu/ml/datasets/default+of+credit+card+clients
    df = pd.read_csv("default of credit card clients.csv", skiprows=1)
    df.columns = [c.lower() for c in df.columns]

    # Drop ID column
    df = df.drop("id", axis=1)
    # change sex and marital status to binary 0/1
    df.sex[df.sex == 2] = 0  # 1 male, 0 female
    df.marriage[df.marriage != 1] = 0
    # column default on next payment is the outcome
    # column marriage is the treatment
    df = df.rename(columns={'default payment next month': 'Y', 'marriage': "A"})

    # rearrange columns so that the outcome and treatment are in columns 1 and 2
    cols = list(df.columns)
    cols.pop(cols.index('A'))
    cols.pop(cols.index('Y'))
    df = df[['Y', 'A'] + cols]

    # create higher level features from the input features
    # theses high level features will be then used in the treatment and outcome functions but will be hidden to the users
    df['age_cycle'] = np.sin(df.age - np.mean(df.age))
    df['risk'] = np.log((df.education + 1) / df.limit_bal)
    f_young = ECDF(df.age)
    df['young'] = f_young(df.age)

    return df


# generate treatment a complex function of measured covariates for observational data
def gen_treatment(df):
    logit = 1 + .2 * df.sex + .6 * df.age_cycle - .25 * np.log(abs(df.bill_amt1) + 1) - 1 * (df.pay_2 < 0)
    p = sigmoid(logit)
    treatments = np.random.binomial(1, p)

    return treatments


# generate outcome as a complex function of measured covariates
# treatment effect is heterogeneous
def gen_outcome(df):
    logit = -1 - .2 * df.risk + .2 * df.age_cycle + 1 * df.sex - .25 * np.log(np.abs(df.bill_amt1) + 1) - 1 * (
            10 ** -5) * df.bill_amt5

    p = sigmoid(df.A * (1 + .2 * df.young + .2 * df.age_cycle) + logit)
    outcome = np.random.binomial(1, p)

    return outcome


def get_experimental_data(df, n, p=0.5):
    df_truncate = df.copy()
    f_age = ECDF(df_truncate.age)
    age_cdf = f_age(df_truncate.age)

    # filter out subpopulation to make the covariate distributions different from observational data
    # remove the first youngest 30% from the population
    df_truncate = df_truncate[age_cdf > 0.3]
    # remove education level > 2 from the population
    df_truncate = df_truncate.query('education <= 2')

    df_exp = df_truncate.sample(n, replace=True)

    df_exp['A'] = np.random.binomial(1, p, n)
    df_exp['Y'] = gen_outcome(df_exp)

    # remove the 'higher level features'
    df_exp = df_exp.drop(['age_cycle', 'risk', 'young'], axis=1)

    return df_exp


def get_observational_data(df, n):
    df_obs = df.sample(n, replace=True)

    df_obs['A'] = gen_treatment(df_obs)
    df_obs['Y'] = gen_outcome(df_obs)
    # remove the 'higher level features'
    df_obs = df_obs.drop(['age_cycle', 'risk', 'young'], axis=1)

    return df_obs


# use this function to generate the ground truth expected outcome for evaluating candidate models
# since the candidate models do not have access to the 'higher level' features, we need to create those features in this function
def get_expected_outcome(df):
    age_cycle = np.sin(df.age - np.mean(df.age))
    risk = np.log((df.education + 1) / df.limit_bal)
    f_young = ECDF(df.age)
    young = f_young(df.age)

    logit = -1 - .2 * risk + .2 * age_cycle + 1 * df.sex - .25 * np.log(np.abs(df.bill_amt1) + 1) - 1 * (
            10 ** -5) * df.bill_amt5

    p = sigmoid(df.A * (1 + .2 * young + .2 * age_cycle) + logit)
    return p


if __name__ == "__main__":
    np.random.seed(1)

    df = gen_population()

    df_exp = get_experimental_data(df, 500, p=0.5)
    df_obs = get_observational_data(df, 10000)

    # save to csv files
    df_exp.to_csv("experimental_data.csv", index=False)
    df_obs.to_csv("observational_data.csv", index=False)
