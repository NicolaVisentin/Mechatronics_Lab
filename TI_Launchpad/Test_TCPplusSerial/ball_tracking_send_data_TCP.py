import cv2
import numpy as np
import socket
import time
import struct

# Set up the TCP/IP server (localhost, port 12345)

server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)   # utilizzo del protocollo TCP
server_socket.bind(('127.0.0.1', 12345))                            # Listen on port 12345
server_socket.listen(1)

print("Waiting for connection...")
client_socket, client_address = server_socket.accept()
print(f"Connection established with {client_address}")

# Dimensioni dell'immagine

frame_width = 640
frame_height = 480

# Crea una maschera statica per selezionare la "zona di interesse" per il rilevamento

mask_static = np.zeros((frame_height, frame_width), dtype=np.uint8)
contour_points = np.array([[0, 70], [335, 70], [335, 420], [0, 420]])   # Coordinate della zona di interesse (rettangolo)
cv2.fillPoly(mask_static, [contour_points], 255)                        # Riempi la zona di interesse con 255 (bianco)

# Acquisizione dalla camera e rilvemento pallina

cap = cv2.VideoCapture(1)   # prova 1 per DroidCam e 0 per camera USB
if not cap.isOpened():
    print("Errore: non riesco ad aprire la videocamera.")
    exit()

while True:
    ret, frame = cap.read()
    if not ret:
        print("Errore: non riesco a leggere il frame.")
        break

    # Converte l'immagine allo spazio colore HSV

    hsv = cv2.cvtColor(frame, cv2.COLOR_BGR2HSV)

    # Definisci l'intervallo per il colore arancione

    lower_orange = np.array([12, 70, 215])   # HSV (H tra 0 e 169, S e V tra 0 e 255)
    upper_orange = np.array([30, 255, 255])  # HSV (H tra 0 e 169, S e V tra 0 e 255)

    # Crea una maschera per il colore arancione

    mask_orange = cv2.inRange(hsv, lower_orange, upper_orange)

    # Applica la maschera statica per selezionare solo la zona di interersse

    mask_combined = cv2.bitwise_and(mask_orange, mask_static)

    # Trova i contorni nella maschera combinata (statica && arancione)

    contours, _ = cv2.findContours(mask_combined, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Tracking pallina

    for contour in contours:
        if cv2.contourArea(contour) > 500:   # Ignore small contours

            # Calculate the bounding rectangle
            x, y, w, h = cv2.boundingRect(contour)

            # Calculate the center coordinates
            center_x = x + w // 2
            center_y = y + h // 2

            # Pack the data as binary (using struct)
            # 'ii' means two integers ('i' is the format character for int in struct)
            message = struct.pack('>HH', center_x, center_y)

            # Send the binary data via TCP socket
            client_socket.sendall(message)

            # Display the coordinates on the video frame
            cv2.putText(frame, f"center_x: {center_x}", (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2)
            cv2.putText(frame, f"center_y: {center_y}", (10, 60), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2)

            # Draw a rectangle around the ball and mark its center
            cv2.rectangle(frame, (x, y), (x + w, y + h), (0, 255, 0), 2)
            cv2.circle(frame, (center_x, center_y), 5, (0, 0, 255), -1)
            
            # Add a sleep to maintain the sampling rate (pause for approximately 1/60 seconds for a 60Hz refresh rate)
            #time.sleep(0.167)

    # Crea una copia del frame originale per aggiungere il contorno bianco alla zona di interesse

    frame_with_contour = frame.copy()
    cv2.polylines(frame_with_contour, [contour_points], isClosed=True, color=(255, 255, 255), thickness=2)

    # Mostra i vari frame nelle finestre

    cv2.imshow("Tracking", frame_with_contour)      # Frame con il tracciamento
    cv2.imshow("Pixel rilevati", mask_combined)     # Maschera combinata

    # Esci premendo 'q'

    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

# Release resources

cap.release()
cv2.destroyAllWindows()
client_socket.close()
server_socket.close()