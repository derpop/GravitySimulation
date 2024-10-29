using System.Collections.Generic;
using UnityEngine;
public interface IObject {
    Vector3 Position { get; set; }
    float mass {get; set;}
}

public class GravityObject : MonoBehaviour, IObject
{
    public float mass {get; set;}       
    public Vector3 velocity = Vector3.zero;   
    public float G;              
    public float softeningFactor; 
    public float dampingFactor;
    public float theta;

    public Vector3 Position { get; set; }

    /// <summary>
    /// Calculates and applies gravitational forces from other particles to update this particle's velocity.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <param name="timeStep">Time step for the Euler integration</param>
    public void ApplyForces(GravityObject[] others, List<PointOctreeNode<GravityObject>> nodes, float timeStep)
    {
        // Calculate total gravitational acceleration acting on this particle
        Vector3 totalForce = CalculateGravitationalForce(others) + CalculateLargeScaleInteractions(nodes);
        Vector3 totalAcceleration = totalForce/mass;
        velocity += totalAcceleration * timeStep;
    }

    /// <summary>
    /// Updates the position of the particle based on its current velocity, using Euler's method.
    /// </summary>
    /// <param name="timeStep">Time step for position update</param>
    public void UpdatePosition(float timeStep)
    {
        // Apply damping to control velocity growth
        // Update position based on the current velocity
        transform.position += velocity * timeStep;
        Position = transform.position;
    }
    public Vector3 CalculateLargeScaleInteractions(List<PointOctreeNode<GravityObject>> nodes){
        Vector3 totalForce = Vector3.zero;

        foreach (var node in nodes) {
            // node.CalculateCenterOfMass();
            if (node.TotalMass == 0) continue; // Skip empty nodes

        Vector3 direction = node.CenterOfMass - this.Position;
        float distance = direction.magnitude;

        if (distance == 0) continue; // Skip if positions overlap to avoid division by zero

        // Only apply the Barnes-Hut criterion without any close-scale calculations
        if ((node.SideLength / distance) < theta) {
            // Treat the node as a single mass point
            float forceMagnitude = G * this.mass * node.TotalMass / (distance * distance);
            totalForce += direction.normalized * forceMagnitude;
        }
    }
    Debug.Log(totalForce);
    return totalForce; // Return the total gravitational force as a Vector3
    }

    /// <summary>
    /// Calculates the total gravitational acceleration on this particle from all other particles.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <returns>Total gravitational acceleration as a Vector3</returns>
    private Vector3 CalculateGravitationalForce(GravityObject[] others)
    {
        Vector3 force = Vector3.zero;

        // Loop through each other particle and calculate the gravitational effect
        foreach (var other in others)
        {
            if (other == this) continue;  // Skip self-interaction

            // Calculate the direction and softened distance between particles
            Vector3 direction = other.transform.position - transform.position;
            float distance = Mathf.Sqrt(direction.sqrMagnitude + softeningFactor * softeningFactor); // Softened distance

            // Calculate gravitational acceleration and add to the total
            force += G * other.mass * this.mass / (distance * distance) * direction.normalized;
        }

        return force;
    }
}
