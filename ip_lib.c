/*
 Created by 31-TUBYSQUAD from 10/04/20 to 30/05/2020.
-Alessandro Giuliani 882483
-Alessandro Celadon 882778
-Biagio Barchielli 881423
-Leonardo De Checchi 882376
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#include<stdlib.h>

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}
/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
/*VERIFIED*/
ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    
    int i,j,k1;
    ip_mat *t;

    /*SPAZIO PER LA STRUCT*/
    t = (ip_mat*)malloc(sizeof(ip_mat));
    /*inizializziamo i valori delle dimesioni della struct*/
    t->h = h;
    t->w = w;
    t->k = k;
    /*inizializziamo stat con il valore v per ogni singolo elemento */ 
    t->stat=(stats*)malloc(sizeof(stats)*k);
    for(k1=0;k1<k;k1++){
        t->stat[k1].max=v;
        t->stat[k1].min=v;
        t->stat[k1].mean=v;
        
    }
    

    /*SPAZIO PER LA PILA DI PUNTATORI*/
    t->data = (float***)malloc(sizeof(float**)*h);

    /*SPAZIO PER LA MATRICE TRIDIMENSIONALE*/
    for(i = 0; i < h; i++){
        t->data[i] = (float**)malloc(sizeof(float*)*w);
        for(j = 0; j < w; j++){
            t->data[i][j] = (float*)malloc(sizeof(float)*k);
        }
    }

    /*INIZIALIZZIAMO LA MATRICE CON IL VALORE V*/
    for(i=0;i<h;i++)
        for(j=0;j<w;j++)
            for(k1=0;k1<k;k1++)
               set_val(t,i,j,k1,v);


    return t;
}

/* Libera la memoria (data, stat e la struttura) */
/*VERIFIED*/
void ip_mat_free(ip_mat *a){
    int i; int j;
    if(a){
    for(i=0;i<a->h;i++){
        for(j=0;j<a->w;j++){
            free(a->data[i][j]); /* andiamo a cancellare la memoria dinamica del puntatore che punta l'array che comprende i 3 valori rgb*/
        }
        free(a->data[i]);/*andiamo a cancellare la memoria dinamica del puntatore doppio ,che punta a sua volta all'array rgb, di cui abbiamo appena fatto la free */ 
    }
    free(a->stat);
    free(a->data);
    free(a);
    }

}

/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats_single(ip_mat *t, int k)
{
    int i, x;
    float cont = 0.0;
    float acc = 0.0, media = 0.0, max, min;

    min = t->data[0][0][k]; /*inizializzo minimo con il k adeguato*/
    max = t->data[0][0][k]; /*inizializzo massimo con il k adeguato*/

    for(i = 0; i < t->h; i++){
        for(x = 0; x < t->w; x++){    /*scorro data in modo da comparare min max e cambiarli se necessario*/
            if (max < t->data[i][x][k])
                max = t->data[i][x][k]; 
            if (min > t-> data[i][x][k])
                min = t->data[i][x][k];
            acc += t-> data[i][x][k]; /*sommo ad acc ogni elemento di data e aumento il cont*/
            cont++;
        }
    }
    media = acc/cont; /*calcolo la media*/ 
    t->stat[k].max = max;
    t->stat[k].min = min;  /*metto all'interno del mio array stats le statistiche*/
    t->stat[k].mean = media;
}

void compute_stats(ip_mat *t)
{
    int k;
    if(!t){
        printf("Puntatore non valido");
        exit(1);
    }
    for(k = 0; k < t->k; k++){
        compute_stats_single(t,k);
    } 
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
/*VERIFIED*/
void ip_mat_init_random(ip_mat *t, float mean, float var)
{
    int x, y, z;
    if(!t){
        printf("Puntatore non valido");
        exit(1);
    }

    for(x = 0; x < t->h; x++){
        for(y = 0; y < t->w; y++){
            for(z = 0; z < t->k; z++){
                t->data[x][y][z] = get_normal_random(mean,var);
            }
        }
    }
    compute_stats(t);
}

/* Crea una copia di una ip_mat e lo restituisce in output */
/*VERIFIED*/
ip_mat * ip_mat_copy(ip_mat * in)
{
    ip_mat *nuova;
    int  x, y, z;
    if(!in){
        printf("Puntatore non valido");
        exit(1);
    }

    nuova=ip_mat_create(in->h,in->w,in->k,0.0);
    

    for(x = 0; x < nuova->h; x++){
        for(y = 0; y < nuova->w; y++){
            for(z = 0; z < nuova->k; z++){
                set_val(nuova, x, y, z, get_val(in, x, y, z));
            }
        }
    }
    compute_stats(nuova);

    return nuova;
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
/*VERIFIED*/
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)
{
    ip_mat *nuova;
    int x, y, z;
    int row=row_end-row_start;
    int col=col_end-col_start;
    if(!t){
        printf("Puntatore non valido");
        exit(1);
    }
    nuova = ip_mat_create(row,col,t->k, 0.0);
    
    /* 2 - 8 row_*/     
    for(x = 0; x < row; x++){
        for(y = 0; y < col; y++){
            for(z = 0; z < t->k; z++){
                set_val(nuova, x, y, z, get_val(t, x+row_start, y+col_start, z));
            }
        }
    }
    compute_stats(nuova);

    return nuova;
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
/*VERIFIED*/
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c)
{
    int i,j,k;

    ip_mat *nuova;
    if(!a){
        printf("Puntatore non valido");
        exit(1);
    }
    nuova = ip_mat_copy(a);
    
    for(i=0 ; i < a->h ; i++)
    {
        for(j=0; j < a->w ; j++)
        {

            for(k=0 ; k < a->k ; k++)
            {

                nuova->data[i][j][k] *= c;

            }
        }
    }
    compute_stats(nuova);

    return nuova;

}
/*VERIFIED*/
/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat * ip_mat_add_scalar(ip_mat *a, float c)
{
    int i,j,k;
    ip_mat *nuova;
    if(!a){
        printf("Puntatore non valido");
        exit(1);
    }
    nuova=ip_mat_copy(a);
    
    for(i=0 ; i < a->h ; i++)
    {
        for(j=0; j < a->w ; j++)
        {

            for(k=0 ; k < a->k ; k++)
            {

                nuova->data[i][j][k] += c;

            }
        }
    }
    compute_stats(nuova);

    return nuova;

}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat *a, ip_mat *b)
{
    if(a->h!=b->h||a->w!=b->w||!a||!b){
        exit(1);
    }
    
    else{
    
        int i,j,k;

        ip_mat *nuova;
        nuova=ip_mat_create(a->h,a->w,a->k,0.0);
        
        for(i=0 ; i < a->h ; i++)
        {
            for(j=0; j < a->w ; j++)
            {

                for(k=0 ; k < a->k ; k++)
                {

                    nuova->data[i][j][k] = (a->data[i][j][k]+b->data[i][j][k])/2.0;
                    
                }
            }
        }
        
        compute_stats(nuova);
        return nuova;
    }
}

/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
 * */
ip_mat *mission_to_the_loot(ip_mat *a,ip_mat *b, int dimensione,ip_mat *t)
{
    
    if(dimensione==0){
         if (a->w==b->w && a->k==b->k)
            return t=ip_mat_create(a->h+b->h,a->w,a->k,0.0);
    }
    else if(dimensione==1){
               
            if (a->h==b->h && a->k==b->k)
                return t=ip_mat_create(a->h,a->w+b->w,a->k,0.0);
    }
    else if(dimensione==2){
               
            if (a->h==b->h && a->w==b->w)
                return t=ip_mat_create(a->h,a->w,a->k+b->k,0.0);   
    }
    else {
            return t;
    }
    printf("Le immagini devono avere almeno due dimensioni concordi oppure almeno uno dei due puntatori non è valido");
    exit(1);
}


ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione)
{
     
    ip_mat *t=NULL;
    t=mission_to_the_loot(a,b,dimensione,t);
    if (t)
    {
        int heigth=0,width=0,dimension=0, i,j,z;
        if (!dimensione)
            heigth=a->h;
        else if (dimensione==1)
                width=a->w;
            else
                dimension=a->k;       
        for(i=0;i<a->h;i++)
            for(j=0;j<a->w;j++)
                for(z=0;z<a->k;z++)
                    {
                        t->data[i][j][z]=a->data[i][j][z];                   
                    }
        for(i=0;i<b->h;i++)
            for(j=0;j<b->w;j++)
                for(z=0;z<b->k;z++)
                    {
                        t->data[i+heigth][j+width][z+dimension]=b->data[i][j][z];                      
                    }
        compute_stats(t);
    }
    return t;
    
}
/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/
/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b)
{   
    if(a->h!=b->h||a->w!=b->w){
        exit(1);
    }
    
    else{
        ip_mat *t;
        int i,j,z;
        t=ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i=0;i<a->h;i++){
            for(j=0;j<a->w;j++){
                for(z=0;z<a->k;z++){
                    t->data[i][j][z]=a->data[i][j][z]+b->data[i][j][z]; /* sommo i data e metto dentro t->data*/
                }
            }
        }
        compute_stats(t); /*creo le stats di t->data*/
        return t; 
    }
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b)
{
    if(a->h!=b->h||a->w!=b->w){
        exit(1);
    }
    
    else{
        ip_mat *t;
        int i,j,z;
        t=ip_mat_create(a->h, a->w, a->k, 0.0);
        for(i=0;i<a->h;i++){
            for(j=0;j<a->w;j++){
                for(z=0;z<a->k;z++){
                    t->data[i][j][z]=a->data[i][j][z]-b->data[i][j][z];
                }
            }
        }
        compute_stats(t); /*creo le stats di t->data*/
        return t; 
    }
}

/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */
ip_mat * ip_mat_to_gray_scale(ip_mat *in)
{
    int x, y, z, l;

    ip_mat * nuova;
    if(!in){
        printf("Puntatore non valido");
        exit(1);
    }
    nuova = ip_mat_create(in->h, in->w, in->k, 0.0);

    for(x = 0; x < in->h; x++){
        for(y = 0; y < in->w; y++){

            float media = 0.0;
            float acc = 0.0;

            for(z = 0; z < in->k; z++){
                acc += in->data[x][y][z];
            }

            media = acc/in->k;

            for(l = 0; l < in->k; l++){
                set_val(nuova, x, y, l, media);
            }
        }
    }
    compute_stats(nuova);

    return nuova;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
/*VERIFIED*/
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){

    if(a->h!=b->h||a->w!=b->w){
        exit(1);
    }
    
    else{
        ip_mat *nuova;
        int i,j,l;
        nuova=ip_mat_create(a->h,a->w,a->k,0.0);

        for(i=0;i<a->h;i++)
            for(j=0;j<a->w;j++)
                for(l=0;l<a->k;l++) 
                    nuova->data[i][j][l]=alpha*(a->data[i][j][l])+(1.0-alpha)*(b->data[i][j][l]);

        compute_stats(nuova);
        return nuova;
    }
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright)
{
    ip_mat *nuova;  
    if(!a){
        printf("Puntatore non valido");
        exit(1);
    } 
    nuova=ip_mat_add_scalar(a,bright);
    return nuova;  

}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
ip_mat *ip_mat_corrupt(ip_mat * a, float amount){
    ip_mat *t;
    ip_mat *nuova;
    int i,j,k;
    if(!a){
        printf("Puntatore non valido");
        exit(1);
    }
    nuova=ip_mat_copy(a);
    t=ip_mat_create(a->h,a->w,a->k,0.0);
    ip_mat_init_random(t,0,1);/*creo un ip mat randomico per combinarsi all'ipmat della foto*/
    
    for(i=0;i<a->h;i++)
        for(j=0;j<a->w;j++)
            for(k=0;k<a->k;k++){
                nuova->data[i][j][k]+=t->data[i][j][k]*amount;/*formula già scritta*/
            }
    clamp(nuova,0,255);
    compute_stats(nuova);
    ip_mat_free(t);
    return nuova;

}

/**** PARTE 3: CONVOLUZIONE E FILTRI *****/

/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f)
{
    ip_mat * padding; ip_mat * nuova;
    int ax, ay, z, fx, fy, pad_h, pad_w;
    if(!a||!f){
        printf("Puntatore non valido dell'immagine o del filtro");
        exit(1);
    }
    if(f->h%2!=0)
        pad_h = (f->h - 1)/2;
    else
        pad_h = (f->h)/2;

    if(f->w%2!=0)
        pad_w = (f->w - 1)/2;
    else
        pad_w = (f->w)/2;


    padding = ip_mat_padding(a, pad_h, pad_w); /* applico il padding all'immagine a */

    nuova = ip_mat_create(a->h, a->w, a->k, 0.0); /* creo un'immagine grande come a da restituire in out */

    for(ax = 0; ax < (padding->h)-2*(pad_h); ax++){

        for(ay = 0; ay < (padding->w)-2*(pad_w); ay++){

            for(z = 0; z < a->k; z++){
                
                float media = 0.0;

                for(fx = 0; fx < f->h; fx++){

                    for(fy = 0; fy < f->w; fy++){

                        media += padding->data[fx+ax][fy+ay][z] * f->data[fx][fy][z];
                    }
                }
                nuova->data[ax][ay][z] = media;
            }
        }
    }

    compute_stats(nuova);
    ip_mat_free(padding);

    return nuova;
}

/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro.
 *
 * I parametri della funzione non subiscono modiche, il risultato viene salvato e restituito in output
 * all'interno di una nuova ip_mat.
 */
ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w)
{
    ip_mat * nuova;
    int x, y, z;
    if(!a){
        printf("Puntatore non valido");
        exit(1);
    }

    nuova = ip_mat_create(a->h + 2*pad_h, a->w + 2*pad_w, a->k, 0.0);

    for(x = 0; x < a->h; x++){
        for(y = 0; y < a->w; y++){
            for(z = 0; z < a->k; z++){
                nuova->data[x+pad_h][y+pad_w][z] = a->data[x][y][z];
            }
        }
    }

    return nuova;
}

/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter(){
    ip_mat *nuova;
    int k;
        nuova=ip_mat_create(3,3,3,0.0);
    for(k=0;k<3;k++){
        nuova->data[0][1][k]=-1.0;
        nuova->data[1][0][k]=-1.0;
        nuova->data[1][1][k]= 5.0;
        nuova->data[1][2][k]=-1.0;
        nuova->data[2][1][k]=-1.0;
    }

    compute_stats(nuova);

    return nuova;
}


/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    ip_mat *nuova;
    int k;
        nuova=ip_mat_create(3,3,3,0.0);
    for(k=0;k<3;k++){
        nuova->data[0][0][k]=-1.0;
        nuova->data[0][1][k]=-1.0;
        nuova->data[0][2][k]=-1.0;
        nuova->data[1][0][k]=-1.0;
        nuova->data[1][1][k]= 8.0;
        nuova->data[1][2][k]=-1.0;
        nuova->data[2][0][k]=-1.0;
        nuova->data[2][1][k]=-1.0;
        nuova->data[2][2][k]=-1.0;
    }

    compute_stats(nuova);

    return nuova;
}

/* Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter(){
    ip_mat*nuova=ip_mat_create(3,3,3,0.0);
    int k;
    for(k=0;k<3;k++){
        nuova->data[0][0][k]=-2.0;
        nuova->data[0][1][k]=-1.0;
        nuova->data[1][0][k]=-1.0;
        nuova->data[1][1][k]= 1.0;
        nuova->data[1][2][k]= 1.0;
        nuova->data[2][1][k]= 1.0;
        nuova->data[2][2][k]= 2.0;
    }

    compute_stats(nuova);

    return nuova;
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(unsigned int w,unsigned int h, unsigned int k){
    float media;
    ip_mat *t;
    media=1.0/(w*h*1.0);
    t=ip_mat_create(w,h,k,media);
    compute_stats(t);
    return t;

}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k1, float sigma){
    int i,j,k,x,y;
    ip_mat *t;
    ip_mat *f;
    float gauss;
    float somma=0;
    int metah=h/2,metaw=w/2;
    t=ip_mat_create(w,h,k1,0.0);
    for(k=0;k<k1;k++){
        somma=0;
        for(i=0;i<w;i++)
            for(j=0;j<w;j++){
                    x=i-metah;
                    y=j-metaw;
                    gauss=(1.0/(2.0*PI*sigma*sigma)*exp(-(x*x+y*y)/(2.0*sigma*sigma)));
                    t->data[i][j][k]=gauss;
                    somma+=gauss;

                }
    }
    f=ip_mat_mul_scalar(t,1.0/somma);
    ip_mat_free(t);
    compute_stats(f);
    return f;
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula valore-min/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat * t, float new_max)
{
    int i,j,z;
    compute_stats(t);

    for(i=0;i<t->h;i++)
    {
        for(j=0;j<t->w;j++)
        {
            for(z=0;z<t->k;z++)
            {
         
                t->data[i][j][z]=(t->data[i][j][z]-t->stat[z].min)*new_max/((t->stat[z].max)-(t->stat[z].min));
        
            }
        }
    }
    compute_stats(t);
}


/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    int i,j,k;
    if(!t){
        printf("Puntatore non valido");
        exit(1);
    }
    for(i=0;i<t->h;i++)
        for (j=0;j<t->w;j++)
            for(k=0;k<t->k;k++)
                if(t->data[i][j][k]<low)
                    t->data[i][j][k]=low;
                else if(t->data[i][j][k]>high)
                    t->data[i][j][k]=high;
    compute_stats(t);
                    
}

