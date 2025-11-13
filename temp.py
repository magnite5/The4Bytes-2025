def analyse(chain_length, molecular_weight, weighted_aromaticity, weighted_instability_index, weighted_pi):
        Strong = 0
        Moderate = 0
        Weak = 0

        analysis_results = [
            analyse_chain_length(chain_length),
            analyse_molecular_weight(molecular_weight),
            analyse_aromaticity(weighted_aromaticity),
            analyse_instability_index(weighted_instability_index),
            analyse_pi(weighted_pi)
        ]

        for value in analysis_results:
            if value == "Strong":
                Strong += 1
            elif value == "Moderate":
                Moderate += 1
            else:
                Weak += 1

        if Strong > Moderate and Strong > Weak:
            return "Strong"
        elif Moderate > Strong and Moderate > Weak:
            return "Moderate"
        elif Weak > Strong and Weak > Moderate:
            return "Weak"
        elif Strong == Moderate and Strong > Weak:
            return "Strong"
        elif Strong == Weak and Strong > Moderate:
            return "Strong"
        elif Strong == Weak and Strong == Moderate:
            return "Strong"
        elif Moderate == Weak and Moderate > Strong:
            return "Moderate"
        else:
            return "Strong"

def analyseValue(value, strong_threshold, moderate_threshold):
    if value > strong_threshold:
        return "Strong"
    elif value >= moderate_threshold:
        return "Moderate"
    else:
        return "Weak"

def analyse_chain_length(chain_length):
    return analyseValue(chain_length, 400, 150)

def analyse_molecular_weight(molecular_weight):
    return analyseValue(molecular_weight, 60000, 20000)

def analyse_aromaticity(weighted_aromaticity):
    return analyseValue(weighted_aromaticity, 0.10, 0.07)

def analyse_instability_index(weighted_instability_index):
    if weighted_instability_index < 30:
        return "Strong"
    elif weighted_instability_index >= 30 and weighted_instability_index <= 40:
        return "Moderate"
    else:
        return "Weak"
    
def analyse_pi(weighted_pi):
    return analyseValue(abs(7 - weighted_pi), 2, 1)
