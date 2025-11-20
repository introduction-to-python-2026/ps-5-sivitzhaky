from string_utils import parse_chemical_reaction, count_atoms_in_reaction
from equation_utils import build_equations, my_solve


def balance_reaction(reaction):
    # 1. parse reaction
    reactants, products = parse_chemical_reaction(reaction)
    reactant_atoms = count_atoms_in_reaction(reactants)
    product_atoms = count_atoms_in_reaction(products)

    # 2. build equation and solve
    equations, coefficients = build_equations(reactant_atoms, product_atoms)
    coefficients = my_solve(equations, coefficients) + [1]

    return coefficients
