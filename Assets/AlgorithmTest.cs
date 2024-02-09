using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class AlgorithmTest : MonoBehaviour
{
    [SerializeField] private int populationSize;
    [SerializeField] private float mutationRate;
    [SerializeField] private float childRatePerGeneration = 1;
    [SerializeField] private float crossoverRate;
    [SerializeField] private int generations;
    [SerializeField] private GameObject lightPrefab; 
    private _Individual bestSolution_Individual;
    private float bestSolutionCost;
    [SerializeField] private List<GameObject> corners;
    private List<_Individual> population = new List<_Individual>();


    //The Problem

    /* The problem for this algorithm is to create the lowest amounts of lights with the best intensity to fully light up 
     * a given room. 
     * 
     * The genetic algorith can be suited very well for this problem as every parmeter for this problem is know and outcomes
     * and changes can be calculated through the aglorithm
     *
     */


    void Start()
    {

        bestSolutionCost = 100000;

        RunGeneticAlgorithm();

        foreach(Chromosome source in bestSolution_Individual.chromosomes)
        {
            if(source.ShouldTurnOn == 1)
            {
                Instantiate(lightPrefab, source.Position, Quaternion.identity);
            }
            
        }
    }

    //Running the algorithm 

    /* This method runs the genetic algorithm by:
     * creating the number of children per generation and creating a list of _Individual of the population size
     * If the new individuals cost is lower than the set starting best solution cost of 100,000 that individual will be set
     * as the best solution and the best soltuion cost will now be the new individuals cost.
     * For each generation children are created from different parents and cross over and mutate and then added to
     * the children list.
     * The children are then added to the population, sorted based on cost and then the best solution is chosen again 
     * based on the new population.
     */

    void RunGeneticAlgorithm()
    {

        int numberOfChildrenPerGen = (int)(childRatePerGeneration * populationSize);

        for(int i = 0; i < populationSize; i++)
        {
            _Individual new_Individual = new _Individual(corners);

            population.Add(new_Individual);

            if(new_Individual.cost < bestSolutionCost)
            {
                bestSolution_Individual = new_Individual;
                bestSolutionCost = new_Individual.cost;
            }
        }

        for(int i = 0; i < generations; i++)
        {
            List<_Individual> children = new List<_Individual>();

            while(children.Count <  numberOfChildrenPerGen) 
            {
                List<_Individual> parents = ChooseParents(population.Count);

                _Individual parent1 = parents[0];
                _Individual parent2 = parents[1];

                List<_Individual> parents_children = parent1.Crossover(parent2,crossoverRate);

                _Individual child1 = parents_children[0];
                _Individual child2 = parents_children[1];

                child1.Mutate(mutationRate);
                child2.Mutate(mutationRate);

                children.Add(child1);
                children.Add(child2);

            }



            foreach(_Individual child in children)
            {
                population.Add(child);
            }


            population.Sort((x, y) => x.cost.CompareTo(y.cost));

            if (population[0].cost < bestSolution_Individual.cost)
            {
                
                bestSolution_Individual = population[0];
                bestSolutionCost = population[0].cost;
                
            }
            Debug.Log(bestSolutionCost);

        }
        Debug.Log(population.Count);
    }

    /* This is a helper method to choose the parents of children and ensures that each parent 
     * is unique for a set of children  
     */

    public List<_Individual> ChooseParents(int numberOfParents)
    {
        _Individual _parent1 = population[UnityEngine.Random.Range(0, numberOfParents)];
        _Individual _parent2 = population[UnityEngine.Random.Range(0, numberOfParents)];

        if (_parent1 == _parent2)
        {
            ChooseParents(numberOfParents);
        }

        List<_Individual> parents = new List<_Individual>();

        parents.Add(_parent1);
        parents.Add(_parent2);

        return parents;
    }
}

//The Individual

/* This is the Individual which has 10 chromosomes */

public class _Individual
{
    internal int maxChromosomes = 10;
    internal List<Chromosome> chromosomes;
    internal float cost;
    internal List<GameObject> corners;
    public _Individual(List<GameObject> _corners)
    {
        corners = _corners;
        chromosomes = new List<Chromosome>();
        for(int i = 0; i < maxChromosomes; i++)
        {
            chromosomes.Add(new Chromosome());
        }

        CalculateCost();
    }

    //Cost Function

    /* The cost for the _Individual is calculated through the the intensity of each light for every corner 
     * The cost is also calculated based on if the chromosone has its ShouldTurnOn set to 1 
     * The cost is calculated by added the cost of each corner from the helper and then adding 1 for everylight turned on
     */

    public void CalculateCost()
    {
        foreach(GameObject corner in corners)
        {
            cost += CostPerCorner(corner);
        }

        //for each light add 1 to cost if shouldturn on is correct 
        foreach(Chromosome chromosome in chromosomes)
        {
            if(chromosome.ShouldTurnOn == 1)
            {
                cost += 1;
            }
        }

    }

    /* This is a helper method to calculate the cost for each corner of the room
     * A ray is sent from each chromosome to the corner to see if the ray hits the light
     * If the ray hits an intensity is calculated based on the formula for the intensity of the light
     * Formula: 1/d^2
     * If the intensity is greater than the threshold which is set as 0.15 1000 is added to the cost for that corner
     */

    public float CostPerCorner(GameObject corner)
    {
        const float intensityThreshold = 0.15f;
        float intensity = 0;
        float cost = 0;
        foreach(Chromosome chromosome in chromosomes)
        {
            RaycastHit hit;
            Ray lightRay = new Ray(chromosome.Position,  (corner.transform.position - chromosome.Position).normalized);

            if(Physics.Raycast(lightRay, out hit))
            {
                if(hit.collider.gameObject == corner.gameObject)
                {
                    intensity+= 1/Mathf.Pow(Vector3.Distance(chromosome.Position, corner.transform.position), 2);
                }
            }

        }

        if(intensity < intensityThreshold)
        {
            cost +=1000;
        }

        return cost;
    }

    /* For mutation of the individual if the mutationRate is greater than a random number under normal distribution 
     * then a choice of the chromosone part to mutate is chosen. Either the position or the byte of if the light should turn on
     * If the mutation is for the position the position of the chromosone is mutated on each axis by a random number between -0.3 and 0.3
     * If the mutation is for ShouldTurnOn the byte will be mutated from 0 to 1 or from 1 to 0
     */

    public void Mutate(float mutateRate)
    {

        foreach (Chromosome chromosome in chromosomes)
        {
            if (RandomNormal(0.0f, 1.0f) < mutateRate)
            {
                int mutateChoice = UnityEngine.Random.Range(0, 2);

                switch(mutateChoice)
                {
                    case 0:
                        {
                            chromosome.Position += new Vector3(UnityEngine.Random.Range(-0.3f, 0.3f), UnityEngine.Random.Range(-0.3f, 0.3f), UnityEngine.Random.Range(-0.3f, 0.3f));
                            break;
                        }
                    case 1:
                        {
                            if(chromosome.ShouldTurnOn == 0)
                            {
                                chromosome.ShouldTurnOn = 1;
                            }
                            else
                            {
                                chromosome.ShouldTurnOn = 0;
                            }
                            break;
                        }
                }
            }
        }

    }

    /* CrossOver that is implmented here is where parent 2's chromosones are set for child 1 and parent 1's chromosones are set as child 2's */

    public List<_Individual> Crossover(_Individual parent2, float crossoverRate)
    {
        _Individual _child1 = this;
        _Individual _child2 = parent2;
        if (RandomNormal(0.0f, 1.0f) < crossoverRate)
        {
            _child1.chromosomes = parent2.chromosomes;
            _child2.chromosomes = this.chromosomes;
        }

        List<_Individual> _children = new List<_Individual>();

        _children.Add(_child1);
        _children.Add(_child2);

        return _children;
    }

    /* RandomNormal acts as creating a random number from normal distribution  */
    public static float RandomNormal(float minValue = 0.0f, float maxValue = 1.0f)
    {
        float u, v, S;

        do
        {
            u = 2.0f * UnityEngine.Random.value - 1.0f;
            v = 2.0f * UnityEngine.Random.value - 1.0f;
            S = u * u + v * v;
        }
        while (S >= 1.0f);

        // Standard Normal Distribution
        float std = u * Mathf.Sqrt(-2.0f * Mathf.Log(S) / S);

        // Normal Distribution centered between the min and max value
        // and clamped following the "three-sigma rule"
        float mean = (minValue + maxValue) / 2.0f;
        float sigma = (maxValue - mean) / 3.0f;
        return Mathf.Clamp(std * sigma + mean, minValue, maxValue);
    }


}

/* The Chromosome for the _Individual represnets a light */

public class Chromosome
{
    public Vector3 Position { get; set; }
    public Byte ShouldTurnOn { get; set; }

    public Chromosome()
    {
        Position = new Vector3(UnityEngine.Random.Range(-4.8f, 4.8f), UnityEngine.Random.Range(-4.8f, 4.8f), UnityEngine.Random.Range(0.5f, 9f));
        ShouldTurnOn = (Byte)UnityEngine.Random.Range(0, 2);
    }
}

// Results and Conclusion

/*this genetic algorigm works as a basis for this problem. The problem can be improved by implmenting rotation tho the chromosome 
 * and more calculations taking into consideration the unity based intensity system and range 
 * 
 * The cost has reached 1000 + the number of lights activated from an example value that I have seen of 8006, this does show improvment of the algorithm.
 * Which does give a base repensentation of where lights should be placed but can be improved with an improvement to the cost function
 * Potentially by making improvements to the positioning and rotations with a spot light and not an area light
 * 
 * 
 */