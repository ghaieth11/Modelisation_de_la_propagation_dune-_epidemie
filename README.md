# Modélisation et Simulation d'une Épidémie

## Description du projet

Dans ce projet, nous nous concentrons sur la modélisation et la simulation de la propagation d'une épidémie à travers une population. Nous avons utilisé des équations différentielles pour modéliser les dynamiques de la population infectée, sensible, rétablie, et en traitement.

Plus spécifiquement, le projet couvre les points suivants :

### Modélisation mathématique :

- Utilisation des équations différentielles pour modéliser la dynamique des populations SIR (Susceptible, Infecté, Rétabli) ainsi qu'une composante supplémentaire T (Traitement).
- Analyse des paramètres influençant la propagation de l’épidémie, tels que le taux de transmission, le taux de guérison, le taux de mortalité, et le taux de transition vers les rétablis.

### Résolution numérique :

- Mise en œuvre de différentes méthodes de résolution numérique :
  - Méthode d'Euler explicite
  - Méthode d'Euler implicite
  - Méthode de Runge-Kutta d'ordre 4
  - Résolution directe via la fonction `solve_ivp` de SciPy
- Comparaison des méthodes pour évaluer leur efficacité et leur précision.

### Simulation en MATLAB :

Une simulation complète a été réalisée en MATLAB pour visualiser l'évolution de l'épidémie au fil du temps et comprendre comment les différentes populations (susceptibles, infectés, rétablis, en traitement) évoluent en fonction des paramètres.

### Rapport théorique et analyse des résultats :

Un rapport détaillé est fourni, qui contient :
- La résolution théorique complète du modèle.
- Des explications approfondies sur les méthodes numériques utilisées.
- Une discussion sur les résultats obtenus et leur pertinence dans le contexte de la modélisation épidémiologique.
- Une analyse critique des limites et des perspectives du modèle.

## Structure du projet

- `src/` : Contient le code source de la simulation Python et des différentes méthodes numériques pour la résolution du système d'équations différentielles.
- `matlab_simulation/` : Contient la simulation réalisée en MATLAB permettant de visualiser la propagation de l'épidémie.
- `report/` : Dossier contenant le rapport détaillé du projet avec la modélisation théorique, les résultats et la discussion.
- `README.md` : Ce fichier décrivant le projet et ses objectifs.

## Utilisation

### Simulation Python

Pour exécuter la simulation en Python et voir la comparaison des différentes méthodes de résolution numérique :

1. Clonez le dépôt GitHub :

    ```bash
    git clone https://github.com/username/nom-du-projet.git
    cd nom-du-projet
    ```

2. Assurez-vous d'avoir installé les dépendances nécessaires :

    ```bash
    pip install numpy scipy matplotlib ipywidgets
    ```

3. Lancez le fichier Python principal pour visualiser les résultats des différentes méthodes numériques :

    ```bash
    python simulation_epidemie.py
    ```

### Simulation MATLAB

Le fichier `epidemie_simulation.m` dans le dossier `matlab_simulation/` permet de visualiser la propagation de l'épidémie.

1. Ouvrez MATLAB et exécutez le fichier :

    ```matlab
    run('epidemie_simulation.m')
    ```

Cela affichera une animation représentant l'évolution des populations dans le temps.

## Contenu supplémentaire

- **Résultats graphiques** : Plusieurs graphiques sont générés pour illustrer l'évolution des populations en fonction des méthodes numériques utilisées.
- **Analyse des erreurs** : Une comparaison des erreurs commises par chaque méthode numérique est également réalisée, afin de déterminer la méthode la plus précise et la plus rapide.

## Conclusion

Ce projet met en œuvre différentes techniques de modélisation et de simulation pour comprendre la propagation d'une épidémie et les facteurs qui influencent cette propagation. Les résultats obtenus fournissent des informations utiles sur la dynamique des maladies et les outils numériques à utiliser dans des scénarios épidémiologiques.

---

**ALOUI Ghaieth**  
Étudiant en Mathématiques Appliquées et Modélisation, POLYTECH Nice Sophia - Université Côte d'Azur
