# Optimization functions
from skopt.space import Real, Integer
from skopt import BayesSearchCV
import xgboost as xgb


# =================================================================
# Bayesian Optimization
# =================================================================
def baye_opti(X_train, y_train):
    print("="*50)
    print("Searching for best model parameters ...\nTraining XGBoost on EXACT 4-vector features")
    print("="*50)

    # Define search space
    search_spaces = {
        'n_estimators': Integer(100, 500),
        'max_depth': Integer(3, 10), 
        #'nthread':
        'learning_rate': Real(0.01, 0.3, 'log-uniform'),
        'subsample': Real(0.6, 1.0),
        'colsample_bytree': Real(0.6, 1.0),
        'min_child_weight': Integer(1, 10),
        'gamma': Real(0, 0.5),
        'reg_alpha': Real(0, 2),
        'reg_lambda': Real(0.5, 3)
    }

    # Bayesian optimization
    bayes_search = BayesSearchCV(
        estimator=xgb.XGBClassifier(
            objective='binary:logistic',
            #eval_metric='auc',
            eval_metric=['auc', 'error'],  # Track both AUC and error rate (1 - accuracy)
            enable_categorical=True,
            random_state=42
        ),
        search_spaces=search_spaces,
        n_iter=10, # Number of optimization steps 50
        scoring='roc_auc',
        cv=5, # 5
        random_state=42,
        n_jobs=-1,
        verbose=1
    )

    bayes_search.fit(X_train, y_train)
    model_best = bayes_search.best_estimator_
    best_params = model_best.get_params()

    print("Best parameters:", bayes_search.best_params_)
    print("Best cross-validation AUC: {:.4f}".format(bayes_search.best_score_))
    
    #return bayes_search, model_best
    return best_params

# =================================================================
# Initialize model parameters
# =================================================================
def set_model_params(X_train, y_train):
    """
    Train classifier using EXACT invariant mass and true opening angle.
    This is MORE ACCURATE and SIMPLER than using Delta R approximations.
    """

    # Features: EXACT physics quantities from 4-vectors
    #features = ['m_gg', 'opening_angle', 'cos_theta', 'E_asym']
    #X = pair_df[features]
    #y = pair_df['is_pi0']


    ## Split
    #X_train, X_test, y_train, y_test = train_test_split(
    #    X, y, test_size = 0.2, random_state = 42
    #)

    # Ultra-simple XGBoost - shallow trees, fast training
    model = xgb.XGBClassifier(
        n_estimators = 100,
        max_depth = 3, # Shallow = fast, interpretable
        learning_rate = 0.1,
        objective = 'binary:logistic',
        #eval_metric = 'auc', # parameter is used in machine learning models (like XGBoost, LightGBM, or sklearn-style APIs) to specify that you want to evaluate your model using AUC (Area Under the ROC Curve) metric.
        eval_metric=['auc', 'error'],  # Track both AUC and error rate (1 - accuracy)
        #use_label_encode=False,
        enable_categorical=True,
        random_state = 42
    )

    params = model.get_params()

    #return model, X_train, y_train, X_test, y_test
    return params