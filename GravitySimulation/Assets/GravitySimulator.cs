using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using Random = UnityEngine.Random;

public class GravitySimulator : MonoBehaviour
{
    public int numParticlesPerGalaxy;    
    public float G;              // Gravitational constant, adjust for simulation scale
    public float softeningFactor; // Softening factor to prevent extreme forces at close distances
    public float dampingFactor;
    public int localRadius;
    public float GalaxySize;
    public float centralMass;
    public Vector2 minMaxMassValues;
    public Vector3 simulationBounds;
    public Vector3 Galaxy1Position;
    public Vector3 Galaxy2Position;
    public Vector3 Galaxy1Velocity;
    public Vector3 Galaxy2Velocity;
    public GravityObject gravityObject;
    private PointOctree<GravityObject> gravityObjects;
    private List<GravityObject> gravList;


    public void CreateGalaxy(int numParticle,Vector3 galaxyPosition, Vector3 initalVelocity){


        // Spawn each particle with a slight offset and initial velocity
        for (int i = 0; i < numParticle; i++)
        {
            // Random distance and angle for each particle within the disk
            float distanceFromCenter = Random.Range(0, GalaxySize);
            float angle = Random.Range(0, 2 * Mathf.PI);

            // Position the particle in a circular disk around the galaxy center
            Vector3 position = galaxyPosition + new Vector3(
                distanceFromCenter * Mathf.Cos(angle),
                Random.Range(-GalaxySize / 20, GalaxySize / 20), // Small vertical variation for disk thickness
                distanceFromCenter * Mathf.Sin(angle)
            );

            // Instantiate particle
            GravityObject particle = Instantiate(gravityObject, position, Quaternion.identity);
            particle.mass = Random.Range(1f, 5f); // Assign random mass for variety

            // Calculate initial orbital velocity for disk rotation
            Vector3 toCenter = galaxyPosition - particle.transform.position;
            Vector3 tangentialVelocity = Vector3.Cross(toCenter.normalized, Vector3.up) * Mathf.Sqrt(G * centralMass / toCenter.magnitude);
            particle.velocity = tangentialVelocity;
            particle.G = G;
            particle.softeningFactor = softeningFactor;
            particle.dampingFactor = dampingFactor;
            particle.velocity += initalVelocity;
            gravList.Add(particle);
        }
    }
    void Start()
    {
        gravList = new List<GravityObject>();
        gravityObjects = new PointOctree<GravityObject>(Math.Max(Math.Max(simulationBounds.x,simulationBounds.y),simulationBounds.z), transform.position, 1);
        CreateGalaxy(numParticlesPerGalaxy,Galaxy1Position,Galaxy1Velocity);        
    }

    // Update is called once per frame
    void Update()
{
    float timeStep = Time.fixedDeltaTime;;
    int numChecks =0;
    gravityObjects = new PointOctree<GravityObject>(
        Math.Max(Math.Max(simulationBounds.x, simulationBounds.y), simulationBounds.z),
        transform.position,
        1
    );
    foreach (GravityObject gravityObject in gravList)
    {
        gravityObjects.Add(gravityObject,gravityObject.transform.position);
    }
    foreach (GravityObject obj in gravList)
    {
        GravityObject[] nearbyObjects = gravityObjects.GetNearby(obj.transform.position, localRadius);
        obj.ApplyForces(nearbyObjects.ToArray(),timeStep);
    }
    foreach (GravityObject gravityObject1 in gravList)
    {
        gravityObject1.UpdatePosition(timeStep);
    }
    Debug.Log(numChecks);
}

    void OnDrawGizmos() 
    {
	gravityObjects.DrawAllBounds();
	gravityObjects.DrawAllObjects();
    }
}
