# Tema 2 Algoritmi avansati - Algoritmti genetici

## Detalii

Proiectul este impartit in 2 parti, fiecare parte rezolva alta problema. Structura proiect:

- datele de intrare se afla in data/settings.txt
- outputul se afla in data/evolution.txt

Pentru ambele probleme populatiile for fi generate random la fiecare rulare.

## Rezolvare gasirea maximului unei functii de gradul 2

Consideram functia de gradul 2

> f(x) = a*x^2 + b*x + c

```txt
population_size -> int
domain_of_definition_start domain_of_definition_end -> float
a b c -> float
recombination_probability -> float
mutation_probability -> float
steps -> int
points_of_crossing -> int
```

Alege daca algoritmul pastreaza membrul elitist setand in constructor parametrul *keep_elitist*.

## Rezolvare knapsack

```txt
population_size -> int
knapsack_weight -> int
number_of_objects -> int
object_max_value -> int
recombination_probability -> float
mutation_probability -> float
steps -> int
points_of_crossing -> int
```

Alege daca algoritmul pastreaza membrul elitist setand in constructor parametrul *keep_elitist*.
