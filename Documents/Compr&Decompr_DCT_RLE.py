import numpy as np
from PIL import Image
import time
import sys


start=time.time()


# Fonction pour la DCT 2D sur un bloc 8x8
def dct_2d(bloc):
    N = 8
    dct_bloc = np.zeros((N, N))

    for u in range(N):
        for v in range(N):
            somme = 0

            for x in range(N):
                for y in range(N):

                    somme += bloc[x, y] * np.cos(((2*x+1)*u*np.pi)/(2*N)) * np.cos(((2*y+1)*v*np.pi)/(2*N))

            C_u = np.sqrt(1/N) if u == 0 else np.sqrt(2/N)
            C_v = np.sqrt(1/N) if v == 0 else np.sqrt(2/N)
            dct_bloc[u][v] = round(C_u * C_v * somme)

    return dct_bloc


# Fonction pour l'inverse DCT 2D sur un bloc 8x8
def idct_2d(bloc):
    N = 8
    idct_bloc = np.zeros((N, N))

    for x in range(N):
        for y in range(N):
            somme = 0

            for u in range(N):
                for v in range(N):

                    C_u = np.sqrt(1/N) if u == 0 else np.sqrt(2/N)
                    C_v = np.sqrt(1/N) if v == 0 else np.sqrt(2/N)
                    somme += C_u * C_v * bloc[u][v] * np.cos(((2*x+1)*u*np.pi)/(2*N)) * np.cos(((2*y+1)*v*np.pi)/(2*N))

            idct_bloc[x, y] = round(somme)

    return idct_bloc




# Fonction de parcours en zigzag pour une matrice 8x8
def zigzag(bloc):
    zigzag_ordre = []
    N = 8
    
    for i in range(2*N - 1):
    
        if i % 2 == 0:
            
            for j in range(i + 1):
                if i - j < N and j < N:
                    zigzag_ordre.append(bloc[i - j, j])
                    
        else:
            for j in range(i + 1):
                if j < N and i - j < N:
                    zigzag_ordre.append(bloc[j, i - j])
                    
    return zigzag_ordre


# Inverse du parcours en zigzag pour un bloc 8x8
def inverse_zigzag(zigzag_liste):
    N = 8
    index = 0
    bloc = np.zeros((N, N))
    
    for i in range(2 * N - 1):
        
        if index >= len(zigzag_liste):  # Vérifie que l'on ne dépasse pas la longueur de la liste
            break
        
        if i % 2 == 0:
            for j in range(i + 1):
                if i - j < N and j < N and index < len(zigzag_liste):
                    bloc[i - j, j] = zigzag_liste[index]
                    index += 1
        else:
            for j in range(i + 1):
                if j < N and i - j < N and index < len(zigzag_liste):
                    bloc[j, i - j] = zigzag_liste[index]
                    index += 1
                    
    return bloc




#Compression/Décompression entropique


#Codage Run-Length (RLE)
def rle(bloc_zigzaggé):
    resultat = []
    zero_compteur = 0
    
    for num in bloc_zigzaggé:
        if num == 0:
            zero_compteur += 1
            
        else:
            resultat.append((int(zero_compteur), int(num)))
            zero_compteur = 0
            
    resultat.append((0, 0))  # Fin de bloc
    return resultat


# Fonction pour décompresser RLE
def irle(rle_compressé):
    resultat = []
    
    for (zeros, num) in rle_compressé:
        resultat.extend([0] * zeros)
        
        if num != 0:
            resultat.append(num)
            
    return resultat







#COMPRESSION
#Bloc + DCT + Quantification + Zigzag + RLE

ind=1
taille_rle=0

#Compression d'un canal
def compresser_canal(canal):
    global facteur_compression, taille_rle
    canal_compressé = []

    for i in range(hauteur_blocs):
        
        sys.stdout.flush()
        t = np.round(i/hauteur_blocs*100,1) 
        sys.stdout.write (str("Compression : [") + str(t) + str("%]") + chr(13))
        time.sleep (0.001)

        for j in range(largeur_blocs):

            # Extraire le bloc 8x8
            bloc = canal[i*8:(i+1)*8, j*8:(j+1)*8]

            # Appliquer la DCT manuellement
            dct_bloc = dct_2d(bloc)

            #Quantifier les coeffs
            bloc_quantifié = np.round(dct_bloc / (facteur_compression * matrice_quantification))

            # Zigzag
            zigzaggé = zigzag(bloc_quantifié)

            # RLE
            rle_liste = rle(zigzaggé)

            taille_rle+=len(rle_liste)

            canal_compressé.append(rle_liste)
    
    return canal_compressé




#DECOMPREESSION
#IRLE + Inverse Zigzag+ Déquantification + IDCT

def decompresser_canal(données_compressées,canal_decompressé):
    global facteur_compression, ind
    données_decompressées=[]

    #IRLE + Inverse zigzag
    for i in range(len(données_compressées)):

        #IRLE
        irle_matrice = irle(données_compressées[i])

        #Inverse zigzag
        inv_zigzag = inverse_zigzag(irle_matrice)


        données_decompressées.append(inv_zigzag)


    #IDCT et déquantification
    bloc_indice=0
    for i in range(hauteur_blocs):
        
        sys.stdout.flush()
        t = np.round(i/hauteur_blocs*100,2) 
        sys.stdout.write (str("Décompression : [") + str(t) + str("%]") + chr(13))
        time.sleep (0.001)
        
        for j in range(largeur_blocs):

            bloc = données_decompressées[bloc_indice] * facteur_compression * matrice_quantification
            
            #IDCT
            idct_bloc = idct_2d(bloc)

            canal_decompressé[i*8:(i+1)*8, j*8:(j+1)*8] = np.clip(idct_bloc, 0, 255)
            bloc_indice += 1
    
    print("Fin décompression",ind)
    ind+=1
    return np.round(canal_decompressé)



# Matrice de quantification standard (luminosité) pour JPEG
matrice_quantification = np.array([
    [16, 11, 10, 16, 24, 40, 51, 61],
    [12, 12, 14, 19, 26, 58, 60, 55],
    [14, 13, 16, 24, 40, 57, 69, 56],
    [14, 17, 22, 29, 51, 87, 80, 62],
    [18, 22, 37, 56, 68, 109, 103, 77],
    [24, 35, 55, 64, 81, 104, 113, 92],
    [49, 64, 78, 87, 103, 121, 120, 101],
    [72, 92, 95, 98, 112, 100, 103, 99] ])




# Charger l'image (en RVB)
image = np.array(Image.open('E:\\TIPE Maths COMPRESSION MP\\perroquet.jpg'))
facteur_compression = 1

# Séparer l'image en trois canaux (R, V, B)
R, V, B = image[:, :, 0], image[:, :, 1], image[:, :, 2]

# Diviser l'image en blocs 8x8
hauteur, largeur = image.shape[0], image.shape[1]
hauteur_blocs = hauteur // 8
largeur_blocs = largeur // 8

print("Hauteur :",hauteur,"/ Largeur :",largeur,"\nTaille image :",hauteur*largeur*3,"pixels = octets")

# Initialiser les matrices pour stocker les canaux compressés
R_zeros = np.zeros((hauteur, largeur))
V_zeros = np.zeros((hauteur, largeur))
B_zeros = np.zeros((hauteur, largeur))


# Compresser et décompresser les trois canaux séparément
R_compressé = compresser_canal(R)
R_decompressé = decompresser_canal(R_compressé,R_zeros)

V_compressé = compresser_canal(V)
V_decompressé = decompresser_canal(V_compressé,V_zeros)

B_compressé = compresser_canal(B)
B_decompressé = decompresser_canal(B_compressé,B_zeros)


'''
fichier = open("compressed_image_RLE.txt", "a")
fichier.write(str(R_compressed))
fichier.close()
'''


# Fusionner les trois canaux compressés
image_decompressée = np.stack([R_decompressé, V_decompressé, B_decompressé], axis=2)


# S'assurer que les valeurs sont dans l'intervalle [0, 255]
image_decompressée = np.clip(image_decompressée, 0, 255).astype(np.uint8)

print("\nTemps d'exécution :",round(time.time() - start,2),"sec")
print("Taille RLE :",taille_rle*2,"octets")

# Sauvegarder l'image compressée
img_decompressée = Image.fromarray(image_decompressée)
img_decompressée.save('decompressed_image_perroquet.jpg')
img_decompressée.show()
