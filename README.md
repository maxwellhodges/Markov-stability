# Markov-stability

for each of the folders Combinatorial and Normalised, cd into the folder and run

python3 setup.py install

the file you want is then in:

build/lib.linux-x86_64-3.5

it will be called something like:

louvainLNL.cpython-35m-x86_64-linux-gnu.so

Copy this file into the top level (Py3Stability in this case, but anywhere you like as long as it's in the 
same folder as stability.py)

NOTE: if you are using the anaconda version of python, you will need to rename the file to:

louvainLNL.so

Do the same for the other folder and then you should have:

stability.py, louvainLCL.so and louvainLNL.so in the same directory

Open an ipython3 notebook and import stability and you're done.  In the provided folder, there is 
an example adjacency matrix and timesteps array for testing.

e.g. for full stability, copy the adjacency matrix and the timesteps into a ipython session and type:

from stability import calculate_full_stability

Then run:

results = calculate_full_stability(adjacency, timesteps)

results is then a pandas dataframe but feel free to play around with the code if you find a different
structure is better for your needs

