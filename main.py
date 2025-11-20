# main.py

# Import the helper functions from the two modules
from string_utils import (
    parse_chemical_reaction,
    count_atoms_in_reaction,
)

from equation_utils import (
    build_equations,
    my_solve,
)


def balance_reaction(reaction):
    """
    Takes a reaction equation as a string (e.g. "Fe2O3 + H2 -> Fe + H2O")
    and returns a list of coefficients that balance the reaction.

    The last coefficient is fixed to 1, and the other coefficients are
    returned as rational numbers (fractions), e.g.:
        "Fe2O3 + H2 -> Fe + H2O"  ->  [1/3, 1, 2/3, 1]
    """

    # 1. Parse the reaction into reactant and product molecules
    reactants, products = parse_chemical_reaction(reaction)

    # 2. Count atoms in each molecule (reactants and products)
    reactant_atoms = count_atoms_in_reaction(reactants)
    product_atoms = count_atoms_in_reaction(products)

    # 3. Build the system of equations and solve for the coefficients
    equations, coefficients = build_equations(reactant_atoms, product_atoms)
    coefficients = my_solve(equations, coefficients) + [1]

    return coefficients
