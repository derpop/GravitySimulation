using System;
using System.Collections.Generic;
using System.Linq;
using Unity.VisualScripting;
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
    public float theta;
    public int refresh = 1;
    public Vector2 minMaxMassValues;
    public Vector3 simulationBounds;
    public Vector3 Galaxy1Position;
    public Vector3 Galaxy2Position;
    public Vector3 Galaxy1Velocity;
    public Vector3 Galaxy2Velocity;
    public GravityObject gravityObject;
    private PointOctree<GravityObject> gravityObjects;
    private List<PointOctreeNode<GravityObject>> nodes;
    private List<GravityObject> gravList;
    private int frameCount;


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
            particle.mass = Random.Range(minMaxMassValues.x, minMaxMassValues.y); // Assign random mass for variety

            // Calculate initial orbital velocity for disk rotation
            Vector3 toCenter = galaxyPosition - particle.transform.position;
            Vector3 tangentialVelocity = Vector3.Cross(toCenter.normalized, Vector3.up) * Mathf.Sqrt(G * centralMass / toCenter.magnitude);
            particle.velocity = tangentialVelocity;
            particle.G = G;
            particle.softeningFactor = softeningFactor;
            particle.dampingFactor = dampingFactor;
            particle.theta = theta;
            particle.velocity += initalVelocity;
            gravList.Add(particle);
        }
    }
    public bool outOfBounds(GravityObject obj) {
    
    Vector3 pos = obj.transform.position;
    
    return pos.x < -simulationBounds.x || pos.x > simulationBounds.x ||
           pos.y < -simulationBounds.y || pos.y > simulationBounds.y ||
           pos.z < -simulationBounds.z || pos.z > simulationBounds.z;
    }

    public void CreateEqualDistribution(int numParticle){
        for (int i = 0; i < numParticle; i++)
        {
            GravityObject particle = Instantiate(gravityObject, new Vector3(Random.Range(-simulationBounds.x/2,simulationBounds.x/2),
                                                                            Random.Range(-simulationBounds.y/2,simulationBounds.y/2),
                                                                            Random.Range(-simulationBounds.z/2,simulationBounds.z/2)),
                                                                            Quaternion.identity);
            particle.mass = Random.Range(minMaxMassValues.x, minMaxMassValues.y);
            particle.G = G;
            particle.softeningFactor = softeningFactor;
            particle.theta = theta;
            particle.dampingFactor = dampingFactor;
            gravList.Add(particle);

        }
    }
    void Start()
    {
        gravList = new List<GravityObject>();
        gravityObjects = new PointOctree<GravityObject>(Math.Max(Math.Max(simulationBounds.x,simulationBounds.y),simulationBounds.z), transform.position, 1);
        CreateGalaxy(numParticlesPerGalaxy,Galaxy1Position,Galaxy1Velocity);
        CreateGalaxy(numParticlesPerGalaxy,Galaxy2Position,Galaxy2Velocity);
        // CreateEqualDistribution(numParticlesPerGalaxy/10);
        nodes = gravityObjects.GetAllNodes();
    }

    // Update is called once per frame
    void Update()
{
    List<GravityObject> toRemove = new List<GravityObject>();
    frameCount++;
    List<PointOctreeNode<GravityObject>> nodes = gravityObjects.GetAllNodes();
    float timeStep = Time.fixedDeltaTime;
    int numChecks =0;
    if(frameCount%refresh==0){
        gravityObjects.rebuild(gravList);
        nodes = gravityObjects.GetAllNodes();
        gravityObjects.getRoot().CalculateCenterOfMass();
        foreach (GravityObject obj in gravList)
        {
        GravityObject[] nearbyObjects = gravityObjects.GetNearby(obj.transform.position, localRadius);
        obj.ApplyForces(nearbyObjects.ToArray(),nodes,timeStep);
        numChecks += nearbyObjects.Length;
        obj.UpdatePosition(timeStep);
        if(outOfBounds(obj)){
            toRemove.Add(obj);
        }
        }
    }
    gravList.RemoveAll(obj => toRemove.Contains(obj));
    foreach (var obj in toRemove) {
        UnityEngine.Object.Destroy(obj.gameObject); // Assuming each GravityObject has an associated GameObject
    }
    
    // Clear toRemove list after destruction
    toRemove.Clear();
    Debug.Log(numChecks);

}

    void OnDrawGizmos() 
    {
        if(gravityObjects!=null){
        gravityObjects.DrawAllBounds();
	    gravityObjects.DrawAllObjects();
        }
    }
}
