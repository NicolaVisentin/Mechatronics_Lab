import cv2
import numpy as np
import serial
import struct
import csv
import time
import datetime
import os

# Configure serial connection

ser = serial.Serial('COM10', 115200, timeout=1)
print(f"Serial connection opened on port {ser.port}")

# Create output file (for saving data)

script_dir = os.path.dirname(os.path.abspath(__file__))  # get script directory
os.chdir(script_dir)

current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
output_file = f"{current_time}.csv"
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['hardware_time', 'python_time', 'x', 'y', 'x_ref', 'y_ref', 'ux', 'uy'])

# Camera acquisition

cap = cv2.VideoCapture(2)   # choose camera: try 2 for DroidCam or 0 for USB camera
if not cap.isOpened():
    print("Error: can't open the selected camera.")
    exit()

# Create a static mask that defines the zone of interest in the frame

frame_width = 640
frame_height = 480

mask_static = np.zeros((frame_height, frame_width), dtype=np.uint8)
contour_points = np.array([[33, 17], [33, 467], [612, 467], [612, 17]])   # Coordinates of the zone of interest (rectangle)
cv2.fillPoly(mask_static, [contour_points], 255)                          # Fill the mask with 255 (white)

# Define HSV interval to detect the ball

lower_HSV = np.array([12, 130, 150])   # HSV (H between 0 ans 169, S and V between 0 and 255)
upper_HSV = np.array([30, 255, 255])   # HSV (H between 0 ans 169, S and V between 0 and 255)

# Acquire and send data at each timestep

start_time = time.time()
while True:

    ret, frame = cap.read()
    if not ret:
        print("Error: can't read the frame.")
        break

    # Get elapsed time
    time_py = time.time() - start_time

    # Convert image in HSV color space
    hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    # Create a mask for the selected HSV bounds and merge it with the static mask
    mask_HSV = cv2.inRange(hsv, lower_HSV, upper_HSV)
    mask_combined = cv2.bitwise_and(mask_HSV, mask_static)

    # Find contours of detected zones (the ones within the static mask and satisfying the HSV mask)
    contours, _ = cv2.findContours(mask_combined, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Tracking and data exchange
    for contour in contours:

        if cv2.contourArea(contour) > 500:    # ignore detected areas that are too small

            ### TRACKING ###

            # Ball position            
            x, y, w, h = cv2.boundingRect(contour)
            center_x = x + w // 2                             # in pixel
            center_y = y + h // 2                             # in pixel
            center_x_mm=round((450)/(612-33)*(center_x-323))  # in mm (NEEDS CALIBRATION)
            center_y_mm=round((350)/(467-17)*(center_y-242))  # in mm (NEDES CALIBRATION)

            # Print data and tracking in the frame
            cv2.putText(frame, f"x: {center_x_mm} mm", (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2)
            cv2.putText(frame, f"y: {center_y_mm} mm", (10, 60), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2)
            cv2.rectangle(frame, (x, y), (x + w, y + h), (0, 255, 0), 2)
            cv2.circle(frame, (center_x, center_y), 5, (0, 0, 255), -1)

            ### SENDING DATA ###

            # Reference position

            x_ref=0
            y_ref=0

            # Prepare binary data package to send:
                # <,> --> litte-endian, big-endian                                                                    
                # f,d --> float single (4 byte = 32 bit), float double (8 byte, 64 bit)
                # b,B --> int8, uint8 (1 byte, 8 bit )
                # h,H --> int16, uint16 (2 byte, 16 bit)
                # i,I --> int32, uint32 (4 byte, 32 bit)
            message = struct.pack('<hhhh', center_x_mm, center_y_mm, x_ref, y_ref) 

            # Send data package via serial                                                        
            ser.write(message)

            ### RECEIVING DATA ###

            try:

                # Wait to receive 12 byte (4 byte for 3 single floats)
                response = ser.read(12)
                if len(response) == 12:

                    # Unpack received data
                    received_data = struct.unpack('<fff', response)  # 3 float (single precision) in little-endian
                    time_hardware, ux, uy = received_data
                    
                    # Print received data in the frame
                    cv2.putText(frame, f"alpha x: {ux:.0f}", (350, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
                    cv2.putText(frame, f"alpha y: {uy:.0f}", (350, 60), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)

            except Exception as e:
                print(f"Error in reading data: {e}")

            ### SAVING DATA ###
            with open(output_file, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([time_hardware, time_py, center_x_mm, center_y_mm, x_ref, y_ref, ux, uy])

    # Plot frame + zone of interest + center of the zone of interest
    cv2.circle(frame, (323, 242), 5, (255, 0, 0), -1)     # center of the zone of interest
    frame_with_contour = frame.copy()
    cv2.polylines(frame_with_contour, [contour_points], isClosed=True, color=(255, 255, 255), thickness=2)

    # Show frames
    cv2.imshow("Tracking", frame_with_contour)
    cv2.imshow("Detected pixels", mask_combined)

    # Quit running code by pressing 'q'
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

# Release resources

cap.release()
cv2.destroyAllWindows()
ser.close()
print("Serial connection closed.")