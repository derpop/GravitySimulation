using System.Collections.Generic;
using UnityEngine;

public class PointOctree<T> where T : IObject{
	public int Count { get; private set; }
	PointOctreeNode<T> rootNode;
	readonly float initialSize;
	readonly float minSize;

	public PointOctree(float initialWorldSize, Vector3 initialWorldPos, float minNodeSize) {
		if (minNodeSize > initialWorldSize) {
			Debug.LogWarning("Minimum node size must be at least as big as the initial world size. Was: " + minNodeSize + " Adjusted to: " + initialWorldSize);
			minNodeSize = initialWorldSize;
		}
		Count = 0;
		initialSize = initialWorldSize;
		minSize = minNodeSize;
		rootNode = new PointOctreeNode<T>(initialSize, minSize, initialWorldPos);
	}
	public void rebuild(List<T> objects){
		rootNode = new PointOctreeNode<T>(initialSize, minSize, Vector3.zero);
		foreach (T obj in objects)
		{
			Add(obj,obj.Position);
		}
	}

	public void Add(T obj, Vector3 objPos) {
		int count = 0;
		while (!rootNode.Add(obj, objPos)) {
			Grow(objPos - rootNode.Center);
			if (++count > 20) {
				Debug.LogError("Aborted Add operation as it seemed to be going on forever (" + (count - 1) + ") attempts at growing the octree.");
				return;
			}
		}
		Count++;
	}

	public bool Remove(T obj) {
		bool removed = rootNode.Remove(obj);
		if (removed) {
			Count--;
			Shrink();
		}
		return removed;
	}

	public bool Remove(T obj, Vector3 objPos) {
		bool removed = rootNode.Remove(obj, objPos);
		if (removed) {
			Count--;
			Shrink();
		}
		return removed;
	}

	public bool GetNearbyNonAlloc(Ray ray, float maxDistance, List<T> nearBy) {
		nearBy.Clear();
		rootNode.GetNearby(ref ray, maxDistance, nearBy);
		return nearBy.Count > 0;
	}

	public T[] GetNearby(Ray ray, float maxDistance) {
		List<T> collidingWith = new List<T>();
		rootNode.GetNearby(ref ray, maxDistance, collidingWith);
		return collidingWith.ToArray();
	}

	public T[] GetNearby(Vector3 position, float maxDistance) {
		List<T> collidingWith = new List<T>();
		rootNode.GetNearby(ref position, maxDistance, collidingWith);
		return collidingWith.ToArray();
	}

	public bool GetNearbyNonAlloc(Vector3 position, float maxDistance, List<T> nearBy) {
		nearBy.Clear();
		rootNode.GetNearby(ref position, maxDistance, nearBy);
		return nearBy.Count > 0;
	}

	public ICollection<T> GetAll() {
		List<T> objects = new List<T>(Count);
		rootNode.GetAll(objects);
		return objects;
	}
	public List<PointOctreeNode<T>> GetAllNodes() {
    List<PointOctreeNode<T>> allNodes = new List<PointOctreeNode<T>>();
    rootNode.GetAllNodesRecursive(allNodes);
    return allNodes;
	}
	public PointOctreeNode<T> getRoot(){
		return rootNode;
	} 



	public void DrawAllBounds() {
		rootNode.DrawAllBounds();
	}

	public void DrawAllObjects() {
		rootNode.DrawAllObjects();
	}

	void Grow(Vector3 direction) {
		int xDirection = direction.x >= 0 ? 1 : -1;
		int yDirection = direction.y >= 0 ? 1 : -1;
		int zDirection = direction.z >= 0 ? 1 : -1;
		PointOctreeNode<T> oldRoot = rootNode;
		float half = rootNode.SideLength / 2;
		float newLength = rootNode.SideLength * 2;
		Vector3 newCenter = rootNode.Center + new Vector3(xDirection * half, yDirection * half, zDirection * half);

		rootNode = new PointOctreeNode<T>(newLength, minSize, newCenter);

		if (oldRoot.HasAnyObjects()) {
			int rootPos = rootNode.BestFitChild(oldRoot.Center);
			PointOctreeNode<T>[] children = new PointOctreeNode<T>[8];
			for (int i = 0; i < 8; i++) {
				if (i == rootPos) {
					children[i] = oldRoot;
				} else {
					xDirection = i % 2 == 0 ? -1 : 1;
					yDirection = i > 3 ? -1 : 1;
					zDirection = (i < 2 || (i > 3 && i < 6)) ? -1 : 1;
					children[i] = new PointOctreeNode<T>(oldRoot.SideLength, minSize, newCenter + new Vector3(xDirection * half, yDirection * half, zDirection * half));
				}
			}
			rootNode.SetChildren(children);
		}
	}

	void Shrink() {
		rootNode = rootNode.ShrinkIfPossible(initialSize);
	}
}