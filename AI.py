import random
import math

# fungsi untuk memecah rumus math
def objectice_function(x1, x2):
    return -(math.sin(x1) * math.cos(x2) + 4/5 * math.exp(1 - math.sqrt(x1**2 + x2**2))) 

# Fungsi untuk membuat kromosom
def create_kromosom():
    return [random.randint(0, 1) for _ in range(8)]

# Fungsi untuk membuat populasi
def create_populasi():
    populasi = []
    for _ in range(8):
        kromosom = create_kromosom()
        populasi.append(kromosom)
    return populasi

def decode_chromosome(chromosome, a, b, n):
    # chromosome adalah kromosom yang ingin didekode
    # a dan b adalah batas bawah dan atas dari interval
    # n adalah jumlah bit yang digunakan untuk merepresentasikan setiap variabel

    x1_bits = chromosome[:n]  # Ambil n bit pertama untuk x1
    x2_bits = chromosome[n:]  # Ambil n bit berikutnya untuk x2

    # Konversi bit ke bilangan desimal
    x1 = a + int("".join(map(str, x1_bits)), 2) * (b - a) / (2 ** n - 1)
    x2 = a + int("".join(map(str, x2_bits)), 2) * (b - a) / (2 ** n - 1)

    return x1, x2

#selection(mencari parent untuk disilangkan)
#dari evaluation tadi dicari yang paling rendah karena kasus min, kalau rendh nilai fit pasti tinggi
def calculate_fitness(objective_values):
    fitness_values = []
    for obj_val in objective_values:
        fitness = 1 / (1 + obj_val) #+1 untuk menghindari pembagian dengan 0
        fitness_values.append(fitness)
    return fitness_values

#hitung nilai probability
def calculate_probability(fitness_values):
    total_fitness = sum(fitness_values)  # Langkah 1: Menghitung total fitness

    probabilities = []
    for fitness in fitness_values:
        probability = fitness / total_fitness  # Langkah 2: Menghitung probabilitas
        probabilities.append(probability)

    return probabilities

#memilih parent dengan roulette whell, harus menghitung cumulative probability values
#Selanjutnya generate nilai random tadi sesuai dengan jumlah kromom
def calculate_cumulative_probability(probabilities):
    cumulative_probabilities = [sum(probabilities[:i+1]) for i in range(len(probabilities))]
    return cumulative_probabilities

def select_chromosomes_using_roulette_wheel(cumulative_probabilities, populasi):
    n = len(populasi)
    new_population = []
    for _ in range(n):
        random_value = random.uniform(0, 1)  # Menghasilkan nilai acak antara 0 dan 1
        selected_chromosome = None
        for i, cumulative_prob in enumerate(cumulative_probabilities):
            if random_value <= cumulative_prob:
                selected_chromosome = populasi[i]
                break
        new_population.append(selected_chromosome)
    return new_population

#crossover dengan rate 25%

def crossover(new_population, crossover_rate):
    random_values = [random.uniform(0, 1) for _ in range(len(new_population))]
    index_kurangdari_rate = [i for i, value in enumerate(random_values) if value < crossover_rate]
    new_population_1 = new_population.copy()  # Create a copy of the original population

    for i in range(0, len(index_kurangdari_rate), 2):
        index1 = index_kurangdari_rate[i]
        if i + 1 < len(index_kurangdari_rate):
            index2 = index_kurangdari_rate[i + 1]
        else:
            # If there's an odd number of selected parents, handle the last one differently
            # You can choose to skip it, for example
            continue
        
        # Perform crossover on selected parents at index1 and index2
        parent1 = new_population[index1]
        parent2 = new_population[index2]

        # Implement crossover logic (you can use a custom crossover function here)
        # For example, let's use a single-point crossover for simplicity
        crossover_point = len(parent1) // 2
        offspring1 = parent1[:crossover_point] + parent2[crossover_point:]
        offspring2 = parent2[:crossover_point] + parent1[crossover_point:]

        # Replace the parents with offspring in the new_population_1
        new_population_1[index1] = offspring1
        new_population_1[index2] = offspring2

    return new_population_1, index_kurangdari_rate

# Mutation function
def mutate_chromosome(chromosome, mutation_rate):
    mutated_chromosome = []
    for bit in chromosome:
        if random.random() < mutation_rate:
            mutated_bit = 1 - bit  # Flip the bit
        else:
            mutated_bit = bit
        mutated_chromosome.append(mutated_bit)
    return mutated_chromosome

# Function to apply mutation to the entire population
def mutate_population(population, mutation_rate):
    mutated_population = []
    for chromosome in population:
        mutated_chromosome = mutate_chromosome(chromosome, mutation_rate)
        mutated_population.append(mutated_chromosome)
    return mutated_population

# MAIN

def main():
    # Inisialisasi variabel
    non_improvement_count = 0
    max_non_improvement = 5  # Jumlah generasi tanpa peningkatan yang diperbolehkan
    generation = 0
    best_fitness_1 = float("-inf")  # Awalnya diatur ke negatif tak terhingga
    best_chromosome_1 = None
    populasi = create_populasi()

    print("+-----------------------------------------------------------------------------------------------+")
    print("| Generation | Best Chromosome                          | Best Fitness      | x1      | x2      |")
    print("+-----------------------------------------------------------------------------------------------+")

    while generation < 100 and non_improvement_count < max_non_improvement:
        a = -10.0  # Batas bawah interval
        b = 10.0  # Batas atas interval
        n = 4  # Jumlah bit yang digunakan untuk merepresentasikan x1 dan x2
        decoded_populasi = []
        nilai_objective = []
        
        for chromosome in populasi:
            x1, x2 = decode_chromosome(chromosome, a, b, n)
            objective_value = objectice_function(x1, x2)
            nilai_objective.append(objective_value)
            decoded_populasi.append((x1, x2))

        # ... Seleksi, crossover, mutasi, dan perhitungan fitness ...
        fitness_values = calculate_fitness(nilai_objective)
        best_fitness = max(fitness_values)
        best_index = fitness_values.index(best_fitness)
        best_chromosome = populasi[best_index]
        probabilities = calculate_probability(fitness_values)
        cumulative_probabilities = calculate_cumulative_probability(probabilities)
        new_population = select_chromosomes_using_roulette_wheel(cumulative_probabilities, populasi)
        crossover_rate = 0.25  # 25%
        new_population_with_crossover,index_kurangdari_rate = crossover(new_population, crossover_rate)
        mutation_rate = 0.25
        new_population_with_crossover = mutate_population(new_population_with_crossover, mutation_rate)

        decoded_populasi = []
        nilai_objective = []
        populasi = new_population_with_crossover 
        for chromosome in populasi:
            x1, x2 = decode_chromosome(chromosome, a, b, n)
            objective_value = objectice_function(x1, x2)
            nilai_objective.append(objective_value)
            decoded_populasi.append((x1, x2))

        fitness_values_new = calculate_fitness(nilai_objective)
        best_fitness_new = max(fitness_values_new)
        best_index_new = fitness_values_new.index(best_fitness_new)
        best_chromosome_new = new_population_with_crossover[best_index_new]
        if best_fitness_new > best_fitness_1:
            best_fitness_1 = best_fitness_new
            best_chromosome_1 = best_chromosome_new
            non_improvement_count = 0  # Reset count
        else:
            non_improvement_count += 1

        x1, x2 = decode_chromosome(best_chromosome_1, a, b, n)
        best_chromosome_1_str = ", ".join(map(str, best_chromosome_1))
        print("| {:<10} | {:<40} | {:<11} | {:<7} | {:<7} |".format(generation, best_chromosome_1_str, best_fitness_1, x1, x2))
        generation += 1

    print("+--------------------------------------+")
    print("Best Chromosome:", best_chromosome_1)
    print("Best Fitness:", best_fitness_1)
    print("x1=", x1)
    print("x2=", x2)

if __name__ == "__main__":
    main()
    
