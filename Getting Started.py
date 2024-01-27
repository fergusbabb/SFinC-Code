from sympy import symbols, Eq, solve, sqrt

# Step 2: Define symbolic variables
x, y, lam, gam = symbols('x y lambda gamma')

# Step 3: Define the dynamical system
f = -3*x+lam*(sqrt(3/2))*y**2 + (3/2)*x*(2*x**2+gam*(1-x**2-y**2))
g = -lam*(sqrt(3/2))*y*x + (3/2)*y*(2*x**2+gam*(1-x**2-y**2))

# Step 4: Set up the system of equations
eq1 = Eq(f, 0)
eq2 = Eq(g, 0)

# Step 5: Solve the system of equations
critical_points = solve((eq1, eq2), (x, y, lam, gam))

print("Critical Points:", critical_points)