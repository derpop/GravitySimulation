using System.Collections.Generic;
using UnityEngine;

public class GravityObject : MonoBehaviour
{
    public float mass;                 // Mass of the particle
    public Vector3 velocity = Vector3.zero;   // Current velocity of the particle
    public float G;              // Gravitational constant, adjust for simulation scale
    public float softeningFactor; // Softening factor to prevent extreme forces at close distances
    public float dampingFactor;    // Damping factor to reduce excessive energy over time

    /// <summary>
    /// Calculates and applies gravitational forces from other particles to update this particle's velocity.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <param name="timeStep">Time step for the Euler integration</param>
    public void ApplyForces(GravityObject[] others, float timeStep)
    {
        // Calculate total gravitational acceleration acting on this particle
        Vector3 totalAcceleration = CalculateGravitationalAcceleration(others);

        // Update velocity based on calculated acceleration
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
    }

    /// <summary>
    /// Calculates the total gravitational acceleration on this particle from all other particles.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <returns>Total gravitational acceleration as a Vector3</returns>
    private Vector3 CalculateGravitationalAcceleration(GravityObject[] others)
    {
        Vector3 acceleration = Vector3.zero;

        // Loop through each other particle and calculate the gravitational effect
        foreach (var other in others)
        {
            if (other == this) continue;  // Skip self-interaction

            // Calculate the direction and softened distance between particles
            Vector3 direction = other.transform.position - transform.position;
            float distance = Mathf.Sqrt(direction.sqrMagnitude + softeningFactor * softeningFactor); // Softened distance

            // Calculate gravitational acceleration and add to the total
            acceleration += G * other.mass / (distance * distance) * direction.normalized;
        }

        return acceleration;
    }
}
