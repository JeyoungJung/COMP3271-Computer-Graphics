#include "MainFrame.h"
#include <iostream>

namespace
{

    float scale = 1.f;
    float aspect = 1.f;

#ifdef __APPLE__
    unsigned int SCR_WIDTH = 600;
    unsigned int SCR_HEIGHT = 600;
#else
    unsigned int SCR_WIDTH = 1000;
    unsigned int SCR_HEIGHT = 1000;
#endif

    void ScrollCallback(GLFWwindow *window, double xoffset, double yoffset)
    {
        scale *= std::pow(1.1f, (float)yoffset);
    }

    void FrameBufferSizeCallback(GLFWwindow *window, int width, int height)
    {
        SCR_WIDTH = width;
        SCR_HEIGHT = height;
        glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);       // Set the viewport to cover the new window
        aspect = (float)SCR_WIDTH / (float)SCR_HEIGHT; // Set the aspect ratio of the clipping area to match the viewport
    }

}

void MainFrame::LeftMouseMove(float start_x, float start_y, float curr_x, float curr_y)
{
    if (modeling_state_ == OBJ_ROTATION)
    {
        // ---------------------------------- Object Rotation ---------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "transform_mat" to rotate the object about its center.

        // Create 2D vectors for start and current mouse positions
        glm::vec2 startPos(start_x, start_y);
        glm::vec2 currPos(curr_x, curr_y);

        // Compute vector from start to current position
        glm::vec2 displacement = currPos - startPos;

        // Compute perpendicular vector to displacement
        glm::vec2 perpDisplacement(-displacement.y, displacement.x);

        // Convert these screen-space vectors into world-space coordinates
        // Subtracting world space position at starting point from world space position at current point gives rotation axis.
        glm::vec3 rotationAxis = glm::normalize(Screen2World(perpDisplacement + startPos) - Screen2World(startPos));

        // Create a rotation matrix that represents small amount of rotation around this axis based on how far mouse has moved.
        float angle = 0.010f * glm::length(perpDisplacement);
        glm::mat4x4 rotationMatrix = glm::rotate(glm::mat4x4(1.f), angle, rotationAxis);

        // Apply this transformation to mesh_
        mesh_.ApplyTransform(rotationMatrix);
    }

    else if (modeling_state_ == OBJ_TRANSLATION)
    {
        // ---------------------------------- Object Translation ------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "trans_mat" to translate the object along the view plane.

        // Convert start and current mouse positions from screen space to world space
        glm::vec3 startPosWorld = Screen2World(start_x, start_y);
        glm::vec3 currPosWorld = Screen2World(curr_x, curr_y);

        // Compute translation vector from start position to current position
        glm::vec3 translationVector = currPosWorld - startPosWorld;

        // Create a transformation matrix that represents translation by this vector.
        glm::mat4x4 translationMatrix = glm::translate(glm::mat4x4(1.f), translationVector);

        // Apply this transformation to mesh_
        mesh_.ApplyTransform(translationMatrix);
    }

    else if (modeling_state_ == OBJ_EXTRUDE)
    {
        // If the current modeling state is set to OBJ_EXTRUDE

        // Define vectors for the origins and directions of two rays:
        // one at the start position and one at the current position
        glm::vec3 rayOrigin0, rayVector0, rayOrigin1, rayVector1, pstart;

        // Compute these rays using Screen2WorldRay function which presumably transforms screen coordinates to world coordinates
        std::tie(rayOrigin0, rayVector0) = Screen2WorldRay(start_x, start_y);
        std::tie(rayOrigin1, rayVector1) = Screen2WorldRay(curr_x, curr_y);

        int face_index = -1;

        // Intersect first ray with mesh faces and return intersected face index and intersection point
        std::tie(face_index, pstart) = mesh_.FaceIntersection(rayOrigin0, rayVector0);

        // Get vertices of intersected face
        const glm::vec3 &vertexA = mesh_.vertices_[mesh_.faces_[face_index][0]];
        const glm::vec3 &vertexB = mesh_.vertices_[mesh_.faces_[face_index][1]];
        const glm::vec3 &vertexC = mesh_.vertices_[mesh_.faces_[face_index][2]];

        // Calculate two edges of the triangle formed by these vertices
        glm::vec3 vec1 = vertexB - vertexA;
        glm::vec3 vec2 = vertexC - vertexA;

        // Calculate normal vector to this triangle (and thus normal vector to this face)
        glm::vec3 normal = glm ::cross(vec1, vec2);

        // Here it seems you're trying to compute an intersection between second (current) mouse position's direction
        // and plane defined by extrusion direction.
        glm ::vec3 M = normal;
        glm ::vec3 R1 = rayVector1;
        glm ::vec3 V = glm ::cross(M, R1);

        // Plane definition. The plane goes through current mouse position in world space
        //(ray origin at current mouse position), and its normal is perpendicular both to M (extrusion direction)
        // and R (current mouse direction).
        glm ::vec3 planePoint = rayOrigin1;
        glm ::vec3 planeNormal = (glm::cross(V, R1));

        // Vector definition. The vector starts from intersection point on initial click,
        // its direction is along extrusion.
        glm ::vec3 vectorStart = (pstart);
        glm ::vec3 vectorDirection = (normal);

        // Intersection computation between defined vector and plane.
        float t = glm::dot(planePoint - vectorStart, planeNormal) / glm :: dot(vectorDirection, planeNormal);

        // Intersection point computation using parameter t.
        glm::vec3 pCurr=vectorStart+t*vectorDirection ;

        // Translation matrix construction that will move initial click point towards computed intersection point,
        //(it will move along extrusion direction).
        glm::mat4x4 transform_mat = glm::translate(glm::mat4x4(1.f),pCurr - pstart);

        // Applying the constructed transformation matrix to the selected face.
        mesh_.ApplyFaceTransform(face_index, transform_mat);
    }
}

void MainFrame::VisualizeWorldSpace()
{

}

// -------------------------------------------------------------------------------------
// -------------------------- No need to change ----------------------------------------
// -------------------------------------------------------------------------------------
void MainFrame::MainLoop()
{
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation, set viewport with width=1000 and height=1000
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "3DModeling", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, FrameBufferSizeCallback);
    glfwSetScrollCallback(window, ScrollCallback);
    // glad: load all OpenGL function pointers
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    const float alpha = 0.3f;
    const float beta = 0.1f;

    const float r = 5.f;
    camera_.LookAt(r * glm::vec3(std::cos(alpha) * std::cos(beta), std::cos(alpha) * std::sin(beta), std::sin(alpha)),
                   glm::vec3(0.f, 0.f, 0.f),
                   glm::vec3(0.f, 0.f, 1.f));

    glEnable(GL_DEPTH_TEST);

    // render loop
    while (!glfwWindowShouldClose(window))
    {
        ProcessInput(window);

        // glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Apply camera projection;
        camera_.Perspective(90.f, aspect, .5f, 10.f);
        camera_.UpdateScale(scale);
        scale = 1.f;
        camera_.ApplyProjection();

        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the display

        DrawScene();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

    // glfw: terminate, clearing addl previously allocated GLFW resources.
    glfwTerminate();
}

void MainFrame::ProcessInput(GLFWwindow *window)
{
    // Key events
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_ROTATION;
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_TRANSLATION;
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_SUBDIVIDE;
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        modeling_state_ = OBJ_EXTRUDE;
    }

    int current_l_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

    // Handle left mouse
    if (current_l_mouse_state == GLFW_PRESS)
    {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (l_mouse_state_ == GLFW_RELEASE)
        {
            LeftMouseClick(xpos, ypos);
            l_click_cursor_x_ = xpos;
            l_click_cursor_y_ = ypos;
        }
        if (l_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f))
        {
            LeftMouseMove(l_click_cursor_x_, l_click_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_l_mouse_state == GLFW_RELEASE)
    {
        if (l_mouse_state_ == GLFW_PRESS)
        {
            LeftMouseRelease();
        }
    }
    l_mouse_state_ = current_l_mouse_state;

    // Handle right mouse
    int current_r_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    if (current_r_mouse_state == GLFW_PRESS)
    {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (r_mouse_state_ == GLFW_RELEASE)
        {
            RightMouseClick(xpos, ypos);
        }
        if (r_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f))
        {
            RightMouseMove(last_cursor_x_, last_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_r_mouse_state == GLFW_RELEASE)
    {
        if (r_mouse_state_ == GLFW_PRESS)
        {
            RightMouseRelease();
        }
    }
    r_mouse_state_ = current_r_mouse_state;
}

void MainFrame::LeftMouseClick(float x, float y)
{
    if (modeling_state_ == OBJ_SUBDIVIDE)
    {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.SubdivideFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
    else if (modeling_state_ == OBJ_EXTRUDE)
    {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.GenExtrudeFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
}

void MainFrame::LeftMouseRelease()
{
    mesh_.CommitTransform();
}

void MainFrame::RightMouseClick(float x, float y)
{
    return;
}

void MainFrame::RightMouseMove(float start_x, float start_y, float curr_x, float curr_y)
{
    glm::vec2 s_start(start_x, start_y);
    glm::vec2 s_cur(curr_x, curr_y);
    glm::vec2 V = s_cur - s_start;
    glm::vec2 A = glm::vec2(-V.y, V.x);
    glm::vec3 rot_axis = glm::normalize(Screen2World(A + s_start) - Screen2World(s_start));
    glm::mat4x4 rot_mat = glm::rotate(glm::mat4x4(1.f), 0.007f * glm::length(A), rot_axis);
    camera_.ApplyTransform(rot_mat);
}

void MainFrame::RightMouseRelease()
{
    return;
}

glm::vec3 MainFrame::Camera2World(const glm::vec3 &x, float w)
{
    return glm::vec3(camera_.view_mat_inv_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::World2Camera(const glm::vec3 &x, float w)
{
    return glm::vec3(camera_.view_mat_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::Screen2World(const glm::vec2 &v, float depth)
{
    float x = v.x / SCR_WIDTH * 2.f - 1.f;
    float y = v.y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * depth;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

glm::vec3 MainFrame::Screen2World(float scr_x, float scr_y, float camera_z)
{
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * -camera_z;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

std::tuple<glm::vec3, glm::vec3> MainFrame::Screen2WorldRay(float scr_x, float scr_y)
{
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec3 o = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 0.f);
    glm::vec3 v = camera_.view_mat_inv_ * v_camera;
    return std::make_tuple(o, v);
}

void MainFrame::DrawScene()
{
    // Draw mesh
    mesh_.Draw();

    VisualizeWorldSpace();
}