using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class PointOctreeNode<T> where T : IObject {
	// Center of this node
	public Vector3 Center { get; private set; }

	// Length of the sides of this node
	public float SideLength { get; private set; }

	// Minimum allowable size for a node
	float minSize;

	// Bounding box of this node
	Bounds bounds = default(Bounds);

	// Objects stored directly in this node
	readonly List<OctreeObject> objects = new List<OctreeObject>();

	// Children nodes, if subdivided
	PointOctreeNode<T>[] children = null;

	// Check if this node has children
	public bool HasChildren { get { return children != null; } }

	// Bounds for each child node
	Bounds[] childBounds;

	// Maximum objects allowed before subdivision
	const int NUM_OBJECTS_ALLOWED = 15;

	// Original bounds size for resizing after temporary expansion
	Vector3 actualBoundsSize;

	public Vector3 CenterOfMass { get; private set; }
	public float TotalMass { get; private set; }


	// Object structure for items in the octree
	class OctreeObject{
		public T Obj;
		public Vector3 Pos;
        public float mass;

    }

	// Constructor
	public PointOctreeNode(float baseLengthVal, float minSizeVal, Vector3 centerVal) {
		SetValues(baseLengthVal, minSizeVal, centerVal);
	}

	// Add an object to this node or its children
	public bool Add(T obj, Vector3 objPos) {
		if (!Encapsulates(bounds, objPos)) {
			return false; // Position not within bounds
		}
		SubAdd(obj);
		return true;
	}

	// Remove an object by reference
	public bool Remove(T obj) {
		bool removed = false;
		for (int i = 0; i < objects.Count; i++) {
			if (objects[i].Obj.Equals(obj)) {
				removed = objects.Remove(objects[i]);
				break;
			}
		}
		if (!removed && children != null) {
			for (int i = 0; i < 8; i++) {
				removed = children[i].Remove(obj);
				if (removed) break;
			}
		}
		if (removed && children != null) {
			if (ShouldMerge()) {
				Merge();
			}
		}
		return removed;
	}

	// Remove an object at a specific position
	public bool Remove(T obj, Vector3 objPos) {
		if (!Encapsulates(bounds, objPos)) {
			return false; // Position not within bounds
		}
		return SubRemove(obj, objPos);
	}

	// Retrieve objects within maxDistance of a ray
	public void GetNearby(ref Ray ray, float maxDistance, List<T> result) {
		bounds.Expand(new Vector3(maxDistance * 2, maxDistance * 2, maxDistance * 2));
		bool intersected = bounds.IntersectRay(ray);
		bounds.size = actualBoundsSize;
		if (!intersected) return;

		// Check contained objects
		for (int i = 0; i < objects.Count; i++) {
			if (SqrDistanceToRay(ray, objects[i].Pos) <= (maxDistance * maxDistance)) {
				result.Add(objects[i].Obj);
			}
		}

		// Check child nodes if subdivided
		if (children != null) {
			for (int i = 0; i < 8; i++) {
				children[i].GetNearby(ref ray, maxDistance, result);
			}
		}
	}

	// Retrieve objects within maxDistance of a position
	public void GetNearby(ref Vector3 position, float maxDistance, List<T> result) {
		float sqrMaxDistance = maxDistance * maxDistance;

		// Check bounds for overlap with search sphere
#if UNITY_2017_1_OR_NEWER
		if ((bounds.ClosestPoint(position) - position).sqrMagnitude > sqrMaxDistance) return;
#else
		bounds.Expand(new Vector3(maxDistance * 2, maxDistance * 2, maxDistance * 2));
		bool contained = bounds.Contains(position);
		bounds.size = actualBoundsSize;
		if (!contained) return;
#endif

		// Check contained objects
		for (int i = 0; i < objects.Count; i++) {
			if ((position - objects[i].Pos).sqrMagnitude <= sqrMaxDistance) {
				result.Add(objects[i].Obj);
			}
		}

		// Check child nodes if subdivided
		if (children != null) {
			for (int i = 0; i < 8; i++) {
				children[i].GetNearby(ref position, maxDistance, result);
			}
		}
	}

	// Get all objects in this node and its children
	public void GetAll(List<T> result) {
		result.AddRange(objects.Select(o => o.Obj));
		if (children != null) {
			for (int i = 0; i < 8; i++) {
				children[i].GetAll(result);
			}
		}
	}

	// Set children for this node
	public void SetChildren(PointOctreeNode<T>[] childOctrees) {
		if (childOctrees.Length != 8) {
			Debug.LogError("Child octree array must be length 8. Was length: " + childOctrees.Length);
			return;
		}
		children = childOctrees;
	}

	// Draw node boundaries
	public void DrawAllBounds(float depth = 0) {
		float tintVal = depth / 7;
		Gizmos.color = new Color(tintVal, 0, 1.0f - tintVal);
		Bounds thisBounds = new Bounds(Center, new Vector3(SideLength, SideLength, SideLength));
		Gizmos.DrawWireCube(thisBounds.center, thisBounds.size);
		if (children != null) {
			depth++;
			for (int i = 0; i < 8; i++) {
				children[i].DrawAllBounds(depth);
			}
		}
		Gizmos.color = Color.white;
	}

	// Draw all objects in this node and children
	public void DrawAllObjects() {
		float tintVal = SideLength / 20;
		Gizmos.color = new Color(0, 1.0f - tintVal, tintVal, 0.25f);
		foreach (OctreeObject obj in objects) {
			Gizmos.DrawIcon(obj.Pos, "marker.tif", true);
		}
		if (children != null) {
			for (int i = 0; i < 8; i++) {
				children[i].DrawAllObjects();
			}
		}
		Gizmos.color = Color.white;
	}

	// Shrink node if possible
	public PointOctreeNode<T> ShrinkIfPossible(float minLength) {
		if (SideLength < (2 * minLength)) return this;
		if (objects.Count == 0 && (children == null || children.Length == 0)) return this;

		int bestFit = -1;
		for (int i = 0; i < objects.Count; i++) {
			OctreeObject curObj = objects[i];
			int newBestFit = BestFitChild(curObj.Pos);
			if (i == 0 || newBestFit == bestFit) {
				if (bestFit < 0) bestFit = newBestFit;
			} else {
				return this;
			}
		}

		if (children != null) {
			bool childHadContent = false;
			for (int i = 0; i < children.Length; i++) {
				if (children[i].HasAnyObjects()) {
					if (childHadContent) return this;
					if (bestFit >= 0 && bestFit != i) return this;
					childHadContent = true;
					bestFit = i;
				}
			}
		}

		return children == null ? this : children[bestFit];
	}

	// Determine which child node would best fit the object
	public int BestFitChild(Vector3 objPos) {
		return (objPos.x <= Center.x ? 0 : 1) + (objPos.y >= Center.y ? 0 : 4) + (objPos.z <= Center.z ? 0 : 2);
	}

    // Check if this node or its children contain any objects
    public bool HasAnyObjects() {
        if (objects.Count > 0) return true;
        if (children != null) {
            for (int i = 0; i < 8; i++) {
                if (children[i].HasAnyObjects()) return true;
            }
        }
        return false;
    }

    // Set up values for the node, including bounds and child bounds
    void SetValues(float baseLengthVal, float minSizeVal, Vector3 centerVal) {
		SideLength = baseLengthVal;
		minSize = minSizeVal;
		Center = centerVal;
		actualBoundsSize = new Vector3(SideLength, SideLength, SideLength);
		bounds = new Bounds(Center, actualBoundsSize);

		float quarter = SideLength / 4f;
		float childActualLength = SideLength / 2;
		Vector3 childActualSize = new Vector3(childActualLength, childActualLength, childActualLength);
		childBounds = new Bounds[8];
		childBounds[0] = new Bounds(Center + new Vector3(-quarter, quarter, -quarter), childActualSize);
		childBounds[1] = new Bounds(Center + new Vector3(quarter, quarter, -quarter), childActualSize);
		childBounds[2] = new Bounds(Center + new Vector3(-quarter, quarter, quarter), childActualSize);
		childBounds[3] = new Bounds(Center + new Vector3(quarter, quarter, quarter), childActualSize);
		childBounds[4] = new Bounds(Center + new Vector3(-quarter, -quarter, -quarter), childActualSize);
		childBounds[5] = new Bounds(Center + new Vector3(quarter, -quarter, -quarter), childActualSize);
		childBounds[6] = new Bounds(Center + new Vector3(-quarter, -quarter, quarter), childActualSize);
		childBounds[7] = new Bounds(Center + new Vector3(quarter, -quarter, quarter), childActualSize);
	}

	void SubAdd(T obj) {
		if (!HasChildren) {
			if (objects.Count < NUM_OBJECTS_ALLOWED || (SideLength / 2) < minSize) {
				OctreeObject newObj = new OctreeObject { Obj = obj, Pos = obj.Position,mass = obj.mass};
				objects.Add(newObj);
				return;
			}
			if (children == null) {
				Split();
				for (int i = objects.Count - 1; i >= 0; i--) {
					OctreeObject existingObj = objects[i];
					int bestFitChild = BestFitChild(existingObj.Pos);
					children[bestFitChild].SubAdd(existingObj.Obj);
					objects.Remove(existingObj);
				}
			}
		}
		int bestFit = BestFitChild(obj.Position);
		children[bestFit].SubAdd(obj);
	}

	bool SubRemove(T obj, Vector3 objPos) {
		bool removed = false;
		for (int i = 0; i < objects.Count; i++) {
			if (objects[i].Obj.Equals(obj)) {
				removed = objects.Remove(objects[i]);
				break;
			}
		}
		if (!removed && children != null) {
			int bestFitChild = BestFitChild(objPos);
			removed = children[bestFitChild].SubRemove(obj, objPos);
		}
		if (removed && children != null) {
			if (ShouldMerge()) {
				Merge();
			}
		}
		return removed;
	}

	void Split() {
		float quarter = SideLength / 4f;
		float newLength = SideLength / 2;
		children = new PointOctreeNode<T>[8];
		children[0] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(-quarter, quarter, -quarter));
		children[1] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(quarter, quarter, -quarter));
		children[2] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(-quarter, quarter, quarter));
		children[3] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(quarter, quarter, quarter));
		children[4] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(-quarter, -quarter, -quarter));
		children[5] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(quarter, -quarter, -quarter));
		children[6] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(-quarter, -quarter, quarter));
		children[7] = new PointOctreeNode<T>(newLength, minSize, Center + new Vector3(quarter, -quarter, quarter));
	}

	void Merge() {
		for (int i = 0; i < 8; i++) {
			PointOctreeNode<T> curChild = children[i];
			for (int j = curChild.objects.Count - 1; j >= 0; j--) {
				objects.Add(curChild.objects[j]);
			}
		}
		children = null;
	}

	static bool Encapsulates(Bounds outerBounds, Vector3 point) {
		return outerBounds.Contains(point);
	}

	bool ShouldMerge() {
		int totalObjects = objects.Count;
		if (children != null) {
			foreach (PointOctreeNode<T> child in children) {
				if (child.children != null) return false;
				totalObjects += child.objects.Count;
			}
		}
		return totalObjects <= NUM_OBJECTS_ALLOWED;
	}
	public void CalculateCenterOfMass() {
    if (!HasChildren) {
        // Leaf node: calculate center of mass from OctreeObject instances directly
        CenterOfMass = Vector3.zero;
        TotalMass = 0f;
        foreach (var obj in objects) {
            CenterOfMass += obj.Pos * obj.mass;  // Use obj.mass directly from OctreeObject
            TotalMass += obj.mass;
        }
        if (TotalMass > 0) CenterOfMass /= TotalMass;
    } else {
        // Non-leaf node: aggregate center of mass and total mass from child nodes
        CenterOfMass = Vector3.zero;
        TotalMass = 0f;
        foreach (var child in children) {
            if (child != null) {
                child.CalculateCenterOfMass();
                CenterOfMass += child.CenterOfMass * child.TotalMass;
                TotalMass += child.TotalMass;
            }
        }
        if (TotalMass > 0) CenterOfMass /= TotalMass;
    }
}

	public static float SqrDistanceToRay(Ray ray, Vector3 point) {
		return Vector3.Cross(ray.direction, point - ray.origin).sqrMagnitude;
	}
	public void GetAllNodesRecursive(List<PointOctreeNode<T>> nodes) {
    // Add this node to the list
    nodes.Add(this);

    // Check if the node has children and recurse if it does
    if (HasChildren) {
        foreach (var child in children) {
            if (child != null) {
                child.GetAllNodesRecursive(nodes);
            }
        }
    }
	}
}
