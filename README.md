# Collection of LP, DP and BGL Problems

## Linear Programming

General tips:

- You can use `LP.set_r` to change the relation of a single constraint. For example `LP.set_r(i, CGAL::LARGER)` will set the relation of constraint i to larger.
- You can turn any LP into a maximization problem by negating the objective function. Don't forget to negate the resulting value again.
- When you have to minize the $L_\infty$ norm of some values $x_1,x_2,...,x_n$, you can use LP. Introduce an unknown $M$ that is bigger than all $|x_1|,|x_2|,...,|x_n|$. I.e. add two constraints $M\geq x_i$ and $M\geq -x_i$ for all $i$. When you then minimize $M$, you effectively minimize the maximum norm.

### Augean Stables

Approach: LP (feasibility check) + Sliding Window

**LP**

We need to check if Hercules can clean the stalls given a certain amount of hours of work on river Alpheius and on Peneius. We don't need to optimize anything, so we don't have to set the objective function. We're only interested in the (in)fesibility, i.e. `sol.is_infeasible()`.

**Sliding Window**

We start with `hours_a = 24` and `hours_p = 0`. If the LP is feasible, i.e. Hercules can clean the stalls, then we can try with less hours on river A, i.e. decrement `hours_a`. Otherwise, we increment `hours_p`. We keep track of the lowest total hours (`hours_a + hours_p`) where the LP is feasible. Furthermore, if the current total hours is greater than the lowest so far, we simply decrement `hours_a` and continue on to the next iteration. If the LP is not feasible at any point, it is impossible for Hercules to clean the stalls.

### The Empire Strikes Back

Approach: Delaunay + LP

Notice that if an asteroid is inside an explosion, the size of that explosion doesn't actually matter. Only when it is outside the radius it will change (to 0). Increasing the radius of an explosion costs nothing. Thus, we first try to maximize each radius without getting caught. We can do this by creating a triangulation of all the bounty hunters and then querying the closest one for each shot. We store the max squared radius for each shot.

Now we can setup an LP. We'll have an unknown for the energy value of each shot. Then, for every asteroid, we sum up the damage each shot contributes. We first check if the asteroid is even within the explosion radius. If not, we skip that shot. Otherwise we compute the inverse of $\max(1, d^2)$ (using K::FT) and set it as the coefficient. For the other side of the inequality, we put the density of the asteroid.

Next, we add the energy constraint. We simply need to sum up all the energy values and set it to be less than the total available energy. Finally, we solve the LP. There is no objective to optimize, we simply want to check if it is feasible.

### Legions

Approach: LP

**Preprocessing**

We need to ensure a consistent orientation of all the lines. We check the sign of the signed distance (we can omit the division by $\sqrt{a^2+b^2}$). If the sign is negative, we negate $a,b,c$. This makes sure that Asterix and Panoramix are on the positive side of each line, i.e. the signed distance will always be positive.

**LP**

We have an unknown non-negative variable $R$ which denotes a kind of "danger radius". The signed distance of our point to all the lines should be greater than $R$, meaning that the point is safe. Furthermore, we can multiply $R$ with the speed $v$ of a legion to incorporate the speed. Instead of the legions moving faster, our danger radius for that legion grows faster. Then
$$
d(p,l)=\frac{ax+by+c}{\sqrt{a^2+b^2}}\geq v\cdot R\\
    ax+by+c\geq v\cdot\sqrt{a^2+b^2}\cdot R\\
    ax+by-v\cdot\sqrt{a^2+b^2}\cdot R\geq -c\\
$$

### Harry Potter

Approach: LP + Binary Search

**LP**

```
Unknowns
    X_ij: (Non-negative) Units going from person i to person j (m variables total)
    S: Max raised suspicion
    
Constraints
    1. Harry has a units => Outflow from Harry must be at most a
    => sum_j X_2j <= a
    
    2. Slughorn must receive at least b units
    => sum_i  e_i1 * X_i1 >=  b
    => sum_i -e_i1 * X_i1 <= -b
    
    3. The total suspicion cannot exceed c
    => sum_ij s_ij * x_ij <= c
    
    4. Sum of outgoing units cannot exceed sum of incoming units, taking into account the efficiency (flow conservation).
        This has to hold for every person (other than Harry and Slughorn)
    => For person i: sum_j e_ji * X_ji >= sum_k X_ik
    =>  sum_k X_ik - sum_j e_ji * X_ji <= 0
    =>   out flow  -     in flow      <= 0
    
    5. S must be greater than the raised suspicion of a transaction ij.
    =>     s_ij * X_ij <= S
    => s_ij * X_ij - S <= 0
    
Constraint IDs
    0: Constraint 1
    1: Constraint 2
    2: Constraint 3
    3 to n: Constraints 4 (aligned with index of person)
    n+1 to n+m: Constraints 5, i.e. n+k+1 for k = 0...m-1
    
Objective
    We want to minimize the max raised suspicion, i.e. minimize S.

To avoid CGAL::Gmpq, we need to rewrite everything with integer values (CGAL:Gmpz)
Lets define Y_ij = X_ij / e''_ij
Then we can replace any X_ij with e''_ij * Y_ij
                    any e_ij * X_ij with e'_ij * Y_ij
```

**Binary Search**

Before anything, we check if the LP is feasible if we use all people. If not, we can skip the entire binary search and optimization.

Otherwise we perform binary search over the number of included people k. To speed up the process, we can only check feasibility at each iteration. Only in the end we optimize the LP corresponding to the lowest feasible k we found. 

### Casterly Rock

Approach: LP

Disclaimer: I've been told that my approach is a bit overcomplicated.

A line (canal) can be modeled with a normal $n$ and X or Y intercept. We can assume that the normal of the sewage canal is $n_s=(1, N)$, since it cannot be horizontal. Since the fresh water canal is orthogonal, its normal is going to be $n_f=(-N, 1)$. We define the X intercept of the sewage canal as $X$, so the sewage canal crosses the X-axis at point $p_s=(X,0)$.

A common house $p_c=(x,y)$ has to lie on the right side of the sewage canal:
$$

n_s\cdot(p_c - p_s)\geq 0\\
(1, N)\cdot(x-X,y)\geq 0\\
x - X + y \cdot N\geq 0\\
- X + y \cdot N\geq x\\
$$

For common houses, we invert the constraint: $X - y \cdot N\geq -x$.

Then, we optimize without an objective to see if these constraints can be satisfied.

Assuming it can be satisfied, we can now compute the length of each sewage pipe, since we know on which side of the canal each house lies. We add up all these lengths and add a constraint $l_1+l_2+...+l_m\leq s$.

To minimize the length of the longest fresh water pipe, we can introduce another unknown $M$ that "captures" this value. We add a constraint $l\leq M$ for each length $l$ and then minimize $M$.


## Dynamic Programming

### Lord Voldemort

Approach: Sliding Window + DP

**Sliding Window**

```
num_breakable_horcruxes[i] = l means that there exists a contiguous interval of l Horcruxes with an evil sum of exactly k that ends at the i-th horcrux.
If no valid interval ends at index i, the value is 0.
    E.g. for k = 3 and evil = [3, 1, 1, 1, 1, 2]
                              |—||———————||————|
                                    |———————|
=> num_breakable_horcruxes = [1, 0, 0, 3, 3, 2]
```

**DP**

```
DP[i][j] stores the maximum number of Horcruxes that can be destroyed with j members while only considering the first i Horcruxes.
E.g. for k = 3 and evil = [3, 1, 1, 1, 1, 2]
                          |—||———————||————|
                                |———————|
- DP[1][1] = 1: One member can destroy the first Horcrux [3].
- DP[4][1] = 3: The member can now destroy the interval [1, 1, 1] instead, which is 3 in total.
- DP[1][2] = 1: Even with a second member, we can still only break the first Horcrux [3], since we only consider that first Horcrux.
- DP[4][2] = 4: One can destroy [3] and the other [1, 1, 1], which is 4 in total.
- DP[6][2] = 5: One can destroy [1, 1, 1] and the other [1, 2], which is 5 in total.

If there is no permissible strategy for i, j, then DP[i][j] should store an invalid state, e.g. -1.
We initialize the DP table with -1.
The base case is DP[i][0] = 0 (0 members can destroy 0 Horcruxes)

Update rules:
1. If j members can break x Horcruxes from the first i-1 Horcruxes,
then they can also break at least x Horcruxes when considering i Horcruxes.

2. Let l = num_breakable_horcruxes[i - 1], then:
if l == 0 (no valid interval ends at i-1):
   There's nothing an additional member can do at this point,
   so it should be at least the same value as with one fewer member, i.e. DP[i][j-1].

if l != 0 (an interval of length l ends at i-1)
   We need to consider the case where the new member breaks the Horcruxes in this interval.
   To not overlap with previous intervals, we jump back l steps, i.e. we look at DP[i-l][j-1].
   Adding the extra Horcruxes we get from the interval, we have DP[i-l][j-1] + l.

Notice that DP[i-l][j-1] + l covers both cases.
```

### Fighting Pits of Meereen

I didn't solve this one. But here's two well documented approaches:

- 6D DP: https://github.com/MariSchn/ETH-Algorithms-Lab-2024/tree/main/Week_13/Fighting_Pits_of_Meereen
- Non-DP: https://github.com/simon-hrabec/Algolab-2020/tree/main/problems/Week%2011%20-%20Fighting%20Pits%20of%20Meereen

### Severus Snape

I didn't solve this one. Here's a well documented approach:

- DP + Greedy: https://github.com/MariSchn/ETH-Algorithms-Lab-2024/tree/main/Week_05/Severus_Snape

### San Francisco

Approach: DP

The game can be modeled as a graph. However, we don't need boost for this problem. Instead, we'll use DP. The DP table is of size n by k + 1, where n is the number of holes and k the max number of moves allowed. DP[u][t] denotes the best score that ends at node u at time t (after t turns). We start at node 0 and time 0, so DP[0][0] = 0 (no score). The other entries are initially invalid, i.e. unreachable.

Next, we'll go turn-by-turn and node-by-node. For a given node, we loop over its outgoing edges and update the entries in the next row (next turn) to be the sum of the current score and the score of the edge. We take the maximum to find the max score. For any leaf node, we can jump back to the start node. This is equivalent to setting the start value to the leaf node's value (still taking the max). If we update DP[u][t+1] and u is a leaf node, then we'll update DP[0][t+1] as well.

If at any point we reach the required score, we can exit early and return the current turn (+1 if we look ahead). If the entire DP table is filled and we didn't reach the required score at any point, then it is impossible to achieve the given score.

We could optimize this by only storing the current and next rows of the DP table, since we never have to do a jump that goes over more than one row. However, it is not needed here, the full DP is already fast enough.

### Burning Coins

Approach: DP

DP[i][j] = Maximum score we can achieve with the first j coins removed and the last i coins removed. If the number of turns (coins) is odd, the starting player (us) gets to pick the last coin. We can fill in those values in the table: DP[n - i - 1][i] = coin_values[i].

Then, for each entry, we compute the scores of taking the first vs. the last coin. If it's our opponents turn, we take the minimum of these scores, since we want to find the max **guaranteed** amount. Otherwise, if it's our turn, we add the value of the picked coin to the repsective scores and take the maximum.

Finally, DP[0][0] holds the maximum amount of money we're guaranteed to win.

### Worm Kingdom

Yeah no clue on this one.

### Alice and the Hurried Rabbit Clan

Approach: DP + DP

**Precomputation (DP)**

To speed up the computation of the second DP, we precompute the total distances of the rows and columns for each entry on the grid. I.e. row_dist[i][j] = total distance of all rabbits in the row of tile i, j (only considering rabbits to the left).

```
Example (single row):
If we have a row with the values (number of rabbits)
      rabbits = [a, b,    c,       d,          e]
Then row_dist = [0, a, 2a+b, 3a+2b+c, 4a+3b+2c+d]

We observe that:
2*a + b = 2a+b
2*(2a+b) + c - a = 3a + 2b + c
2*(3a+2b+c) + d - (2a+b) = 4a + 3b + 2c + d

=> row_dist[i+1] = 2*row_dist[i] + rabbits[i] - row_dist[i-1]

We already know that the first row is all zeros, and the second row is the first row of the rabbits array.
This is analogous for columns.
```

**Main DP**

DP[i][j] := Min distance to position (i, j) (starting from (0, 0)). Since we can only move right or down, DP[i][j] only depends on DP[i-1][j] and DP[i][j-1], because every move will only add a distance related to that move (independent of the previous path).

```
// Option A: Consider left path as previous path

// We need to add the distances of the new rabbits
// These are in column j, above the tile (i, j)
// i.e. rabbits (0, j), (1, j), ... (i-1, j)

int64_t left_path_dist = DP[i][j - 1] + col_dist[i][j];

// Option B: Consider top path as previous path

// We need to add the distances of the new rabbits
// These are in row i, to the left of the tile (i, j)
// i.e. rabbits (i, 0), (i, 1), ... (i, j-1)
int64_t top_path_dist = DP[i - 1][j] + row_dist[i][j];

// Finally, take the better option, i.e. the minimum
DP[i][j] = std::min(left_path_dist, top_path_dist);
```

### DHL

Approach: DP

First, consider decrementing all weights and volumes by 1. This will simplify $(S_a - k_a)\cdot(S_b - k_b)$ into $S_a'\cdot S_b'$.

We observe that it is never beneficial to pick more than 1 parcel from both the stacks. In each turn, we pick exactly 1 parcel from a stack, and some (>= 1) from the other. This simplifies the states we need to consider in each iteration.

DP[i][j] denotes the minimum cost to process the first i parcels from stack A and the first j parcels from stack B. To compute it, we consider the three cases:
1. DP[i-1][j-1]: A new truck arrives. It takes exactly one parcel from A and one from B. 
2. DP[i-1][j]: The current truck stays. It takes an additional parcel from stack A.
3. DP[i][j-1]: The current truck stays. It takes an additional parcel from stack B.

We simply take the best of those options (minimum), and then add $A_i\cdot B_i$ to it (values need to have been decremented by 1 for this!). You can verify that in all three cases the cost increases by exactly this amount.

Finally, the solution is DP[n][n].

## BGL

### Planet Express

Approach: Strongly Connected Components + Dijkstra's

To figure out which planets are linked, we compute the strongly connected components. Next, we compute the number of planets that are part of the teleportation network for each SCC. The cost of teleportation is exactly one less than that (for any two teleportation planets inside that SCC).

To include the teleportation connections, we introduce a new "hub" node for each SCC. We connect all teleportation planets to their corresponding hub (via the SCC map), in both directions. The edge going into the hub has a weight equal to the cost, and the outgoing edge has weight 0.

We add another node that connects to each warehouse with weight 0. We can then perform Dijkstra's from this new node to effectively get the shortest path starting from any of the warehouses. If the shortest path to our target is longer than one million (1 second in microseconds), then we cannot complete the delivery.

### Tiles

Approach: Matching

We parse the input and create a graph of the tiles. Each tile "." gets connected to its four neighbor tiles (if they are to be tiled "."). We then compute a maximum matching of this graph. Each match corresponds to a domino, since it covers exactly two adjacent tiles. Thus, the maximum number of tiles that can be covered is twice the maximum matching size. If it is not equal to the number of tiles ".", then we can't tile the whole garden.

### Placing Knights

Approach: Matching + König Theorem

We construct a graph of all present tiles of the chessboard. We connect vertices which are a "knight's move" apart. We want to find a maximum subset of tiles where no two knights can threaten each other. In other words, we want to find a maximum subset of vertices, where no two vertices are connected. This is exactly the definition of the maximum independent set.

**Theorem (König)**

In a bipartite graph, the number of edges in a maximum matching is equal to the number of vertices in a minimum vertex cover.

When a knight moves in chess, it color of the tile it stands on changes. Therefore, each vertex in our graph is only connected to vertices corresponding to tiles of the other color. This means that our graph is bipartite, so we can use König's theorem to find the number of vertices in a minimum vertex cover.

Finally, remember that the maximum independent set is the complement of the minimum vertex cover: MaxIS = V \ MinVC. Therefore, we can very easily compute its size.

### Targaryen

Approach: Dijkstra's + Matching

We first setup a weighted graph as described in the problem statement. We again have a multi-source shortest path problem (like in Planet Express), so we add another node to all our sources (barracks) without any edge weights. We then perform Dijkstra's to find out the shortest distance to each vertex from any of the barracks.

Now, we construct a new graph consisting only of the vertices that are within reach of a barracks. Making a road safe requires building a barricade at both ends (vertices). We can thus think of a safe road as an edge that covers two vertices. For regular intersections, we can only build a single barricade. Thinking again of edges, this means that two edges corresponding to safe roads cannot cover the same vertex. This is exactly what a matching is.

To incorporate plazas, where we can build two barricades on incident roads, we have to modify the graph a bit. For any plaza node (in range), we connect it normally to its neighbor. We add a duplicate of the plaza node and also connect it to the same neighbor. This ensures that we have two choices for the plaza node, but still only one for the other. Also remember that no two plazas are incident, so we don't have to handle that case.

Finally, to compute the maximum number of safe roads, we compute the maximum matching.

### Return of the Jedi

Approach: Union Find

We're given a description of an MST algorithm and we have to find the smallest spanning tree that is different (by at least one edge) from the MST produced by that algorithm.

We start by computing that MST using a union find datastructure. Then, for each edge in the MST, we construct a new MST without using that edge. Again we can use union find to do this efficiently. The solution is the minimum weight of all of these newly constructed MSTs, since by construction, all of them differ by at least one edge.
