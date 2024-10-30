using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Unity.Burst;
using Unity.Collections;
using System;
public interface IObject {
    Vector3 Position { get; set; }
    float mass {get; set;}
}
public struct GravityObjectData{
    public float3 position;
    public float3 velocity;
    public float mass;
}
public struct NodeData{
    public float3 CenterOfMass;
    public float TotalMass;
    public float3 position;
    public float Size;
    public bool isLeaf;
}
public class GravityObject : MonoBehaviour, IObject
{
    public float mass {get; set;}       
    public Vector3 velocity = Vector3.zero;     
    public Vector3 Position { get; set; }

    /// <summary>
    /// Calculates and applies gravitational forces from other particles to update this particle's velocity.
    /// </summary>
    /// <param name="others">List of other GravityObject particles in the simulation</param>
    /// <param name="timeStep">Time step for the Euler integration</param>
    public void ApplyForces(Vector3 force,float timeStep)
    {
        // Calculate total gravitational acceleration acting on this particle
        Vector3 totalAcceleration = force/mass;
        velocity += totalAcceleration * timeStep;
    }

    /// <summary>
    /// Updates the position of the particle based on its current velocity, using Euler's method.
    /// </summary>
    /// <param name="timeStep">Time step for position update</param>
    public void UpdatePosition(float3 newPos)
    {
        transform.position = newPos;
        Position = transform.position;
    }
    public Vector3 getPosition(){
        Position = transform.position;
        return Position;
    }
}
