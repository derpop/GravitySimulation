using UnityEngine;

public class CameraController : MonoBehaviour
{
    public float movementSpeed = 10f;   // Speed of movement
    public float rotationSpeed = 100f;  // Speed of rotation

    private float yaw = 0.0f;
    private float pitch = 0.0f;

    private void Update()
    {
        HandleMovement();

        // Only rotate the camera if the right mouse button is held down
        if (Input.GetMouseButton(1)) // 1 is the right mouse button
        {
            // Lock the cursor while rotating
            Cursor.lockState = CursorLockMode.Locked;
            HandleRotation();
        }
        else
        {
            // Unlock the cursor when not rotating
            Cursor.lockState = CursorLockMode.None;
        }
    }

    private void HandleMovement()
    {
        // Get input for movement
        float horizontal = Input.GetAxis("Horizontal"); // A/D or Left/Right arrow
        float vertical = Input.GetAxis("Vertical");     // W/S or Up/Down arrow
        float upDown = 0f;

        // Use Q/E for moving up and down
        if (Input.GetKey(KeyCode.Q)) upDown = -1f;
        if (Input.GetKey(KeyCode.E)) upDown = 1f;

        // Calculate movement direction
        Vector3 direction = new Vector3(horizontal, upDown, vertical).normalized;

        // Move the camera
        transform.Translate(direction * movementSpeed * Time.deltaTime);
    }

    private void HandleRotation()
    {
        // Get mouse input for rotation
        yaw += rotationSpeed * Input.GetAxis("Mouse X") * Time.deltaTime;
        pitch -= rotationSpeed * Input.GetAxis("Mouse Y") * Time.deltaTime;

        // Limit pitch to avoid flipping the camera
        pitch = Mathf.Clamp(pitch, -90f, 90f);

        // Apply rotation
        transform.eulerAngles = new Vector3(pitch, yaw, 0.0f);
    }
}
