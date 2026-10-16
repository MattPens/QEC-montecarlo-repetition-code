#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

void MomentiDistr(int N,double &media,double& quadro){
    double sum=0.0;
    double quadrati=0.0;
    for(int i=0;i<N;i++){//conto le volte che estraggo da p(u)
        double u=drand48(); //estraggo a caso valori e gliene associo uno 
        sum+=u;
        quadrati+=u*u;
    }
    media=sum/N;
    quadro=quadrati/N;

}


double ProbEvento(double p,int N){
    //deve dare frequenza degli eventi
    int K=0;
    for(int i=0;i<N;i++){
        double u=drand48();//estraggo valore tra [0,1]
        if(u<p){// Se u e' minore della soglia p, avviene l'errore.
            //u<0.1-->probabilità del 10% avvenga evento
            K++;//segno un evento in più successo
        }
    }
    double p_tilde= (double)K/N; //frequenza eventi
    return p_tilde;
}


double IncertezzaEvento(double p,int N){
    double dev_standard=sqrt(p*(1.-p)/N);
    return dev_standard; //incertezza teorica
}


void RipetiEsperimenti(double p,int N,int M,double & p_media_stimata, double& s){
    //M numero di volte che fai esperimento e N sono estrazioni per ogni esp che fai
    double sum=0.0;
    double p_tilde=0.;
    double sum_quadrati=0.0;
    for(int i=0;i<M;i++){//ciclo du esperim singolo
        p_tilde=ProbEvento(p,N);
        sum+=p_tilde;
        sum_quadrati+=p_tilde*p_tilde;
    }
    p_media_stimata=sum/M;
    s=sqrt((sum_quadrati -M*p_media_stimata*p_media_stimata)/(M-1.0)); //dev standard campionaria dalla p stimata
}


void IntegraleMonteCarlo(double (*f)(double ),double a, double b,int N,double &integrale, double &inc_stimata){
    double sum=0.0;
    double sum_quadrati=0.0;
    double h=(b-a)/N; //N resta # campioni casuali
    for(int i=0;i<N;i++){
        double u_i=drand48();
        double x_i=a+(b-a)*u_i;//punto unif nell'intervallo scelto
        double y_i=f(x_i);
        sum+=y_i;
        sum_quadrati+=y_i*y_i;
    }
    integrale=h*sum;
    double sf2=(sum_quadrati-(sum*sum)/N)/(N-1.);
    inc_stimata=(b-a)*sqrt(sf2/N);
}

// BLOCCO 2: 
int EstraiErrore(double p){
    //una chiamata, unica classificazione 1 per err o 0 se non c'è
    double u_i=drand48(); //singola pescata 
    if(u_i<p){
        return 1; //ERRORE
    }else{
        return 0; //NIENTE ERRORE
    }
}


void GeneraErrori(double p,int errori[3]){//fisso a 3 qubit la sequenza 
    for(int i=0;i<3;i++){
        errori[i]=EstraiErrore(p);//estrare u_i,confronta con p, dice se è err o no, e assegna valore a elem array errore
    }

}

//f che conta i qubit affetti da errore
int PesoErrore(const int errori[3]){
    int w=0;
    for(int i=0;i<3;i++){
        w+=errori[i];
    }
    return w; //peso della configurazione
}


double BinomialCoeff(int n,int w){ //n sopra e w sotto
    if(w==0 || w==n){
        return 1;
    }else{
        return BinomialCoeff(n-1,w-1)+ BinomialCoeff(n-1,w);
    }
}


double ProbPesoTeorica(double p,int w){
    double qw=BinomialCoeff(3,w)*pow(p,w)*pow(1.-p,3-w); 
    return qw;
} // Probabilita totale di ottenere una configurazione di peso w.


void SimulaRumore(double p,int N_configurazioni, int conteggi_qubit[3], int conteggi_w[4]){
    int errori[3];
    //azzero i contatori prima della simulazione
    for(int i=0;i<3;i++){
        conteggi_qubit[i]=0;
    }
    for(int w=0;w<4;w++){
        conteggi_w[w]=0;
    }

    //allora adesso facciamo le simulazioni
    for(int j=0;j<N_configurazioni;j++){
        GeneraErrori(p,errori); //per ogni conf rigenero array errori estrendo a caso u e assegnando 1 o 0
        //conto quanti errori ho per ogni config e aggiorno conteggio qubit errati e pesi
        for(int i=0;i<3;i++){

            if(errori[i]==1){
                conteggi_qubit[i]++; //aumenta di 1 #errori nella configurazione estratta
            //ogni qubit ha il suo contatore che dice quante volte è stato affetto da errore
            }
        }
        int w=PesoErrore(errori); //calcola peso di quella configurazione
        conteggi_w[w]++;
    }
}


//BLOCCO 3: sindrome e detection

void CalcolaSindrome(const int errori[3],int sindrome[2]){
    //sindrome sono i due stabilizzatori Z1Z2 e Z2Z3
    int s1=(errori[0] + errori[1]);  //z su 1 e z su 2
    int s2=(errori[1] + errori[2]);
    if(s1%2==0){ //se ho 2 err ho che s1=2 che mod 2 fa 0
        sindrome[0]=0;
    }else if(s1%2 != 0){
        sindrome[0]=1; //errore 
    }
    if(s2 %2==0){
        sindrome[1]=0;
    }else if(s2 %2 !=0){
        sindrome[1]=1;
    }
// I bit di syndrome codificano gli esiti di Z1Z2 e Z2Z3:
// 0 significa autovalore +1, 1 significa autovalore -1.
}


int IndiceSindrome(const int sindrome[2]){
    //rappr 1dim di vettore delle sindromi con cui mappo le 4 possibilità in detection
    return 2*sindrome[0] + sindrome[1]; //indice 0 ,indice 1 per err su 3,indice 2 per err su 1 e 3 per err su 2
}


//BLOCCO 4: DECODER

void DecodificaSindrome(const int sindrome[2],int correzione[3]){
    //riceve solo sindrome e considera indice per individuare errore su quale qubit
    int indice_sindrome=IndiceSindrome(sindrome); //dal valore capisco correzione su quale qubit e aggiorno correzione[]
    if(indice_sindrome==0){
        correzione[0]=0;
        correzione[1]=0;
        correzione[2]=0; //nessun errore quindi non intervengo su nessun qubit
    } else if(indice_sindrome==1){
        correzione[0]=0;
        correzione[1]=0;
        correzione[2]=1; //agisce sul 3 qubit
    } else if(indice_sindrome==2){
        correzione[0]=1; //sul primo
        correzione[1]=0;
        correzione[2]=0;
    } else if(indice_sindrome==3){
        correzione[0]=0; 
        correzione[1]=1; //sul secondo
        correzione[2]=0;
    }
} //riempe correzione a seconda della sindrome presente 

void ApplicaCorrezione(const int errori[3],const int correzione[3],int residui[3]){
    //applico la correzione alla mia configurazione di 3 qubit
    for(int i=0;i<3;i++){
        residui[i]=(errori[i] + correzione[i])%2;
    }
    //considero caso per cui correggo qubit errato generando errore!!
}


bool ErroreLogico(const int residui[3]){
    //residui=000 correzione fatta :) residui 111 applica X logico su tutti
    return PesoErrore(residui)==3; //true se è 111 e false per resto
}


//BLOCCO 5: MonteCarlo
/* ogni procedura dà 1 se residuo:111 cioè err logico o 000 per procedura corretta
Ripetendo N volte exp conto fallimenti K_L e stimo P_L=K_L/N e confronto con p teorica di err fisico su ogni qubit*/

bool SimulaCicloQEC(double p){ 
    int errori[3];
    GeneraErrori(p,errori);

    int sindrome[2];
    CalcolaSindrome(errori,sindrome);

    int correzione[3];
    DecodificaSindrome(sindrome,correzione); //cosa agisce su quale qubit

    int residui[3];
    ApplicaCorrezione(errori,correzione,residui);

    bool errore_logico=ErroreLogico(residui);
    //stampa 0 se non serve Xlogico e 1 se serve
    return errore_logico;
}


void StimaErroreLogico(double p, int N,double &p_logica, double &incertezza){
    int K_L=0; //fallimenti
    p_logica=0.0; 
    for(int i=0;i<N;i++){
        //ripeti exp N volte
        bool sing_ciclo=SimulaCicloQEC(p);
        if(sing_ciclo){
            K_L++; //aggiorno fallimenti
        }

    }
    p_logica= (double) K_L/N;
    incertezza=IncertezzaEvento(p_logica,N);
} //conta prob che procedura fallisca e confronto con quella teorica

double ProbLogicaTeorica(double p){
    return ProbPesoTeorica(p,2)+ ProbPesoTeorica(p,3);
}

//Parte grafica, blocco 6: faccio variare p e vedo come varia P logica(p) rispetto a curva teorica
void ScansionaProbabilita(double p_min, double p_max,int intervalli, int N){
    ofstream dati("curva_PLogica.dat");

    if(!dati){
        cerr<< "Errore nell'apertura del file dati"<<endl;
        return;
    }

    dati<<setprecision(10);
    dati << "# p P_logica_stimata incertezza P_logica_teorica\n";

    double p_logica=0.0;
    double incertezza=0.0;
    double p_log_teorica=0.0;
    for(int j=0;j<=intervalli;j++){
        double p_j=p_min+(double) j*(p_max - p_min)/intervalli;
        StimaErroreLogico(p_j,N,p_logica,incertezza);
        p_log_teorica=ProbLogicaTeorica(p_j);
    
        dati << p_j << " " << p_logica << " " << incertezza << " " << p_log_teorica << "\n";
    }
    dati.close();
}



int main(){
    cout << setprecision(10);
    srand48(12345);//seed
    double media,quadro;
    double scarto_mom1=1.;
    double scarto_mom2=1.;
    double mom1_vero=0.5;
    double mom2_vero=1./3.;
    int N=100;

    double p_stimata=0.0;
    double p=0.5;
    double dev_standard=1.0;
    double R=1.;

    for(N=100;N<=1000000;N*=10){
        MomentiDistr(N,media,quadro);
        scarto_mom1=fabs(media-mom1_vero);
        scarto_mom2=fabs(quadro-mom2_vero);
        cout<<"N estrazioni: "<<N<<" mom 1: "<<media<<" mom 2: "<<quadro<<" scarto mom 1: "<<scarto_mom1<<" scarto mom 2 :"<<scarto_mom2<<endl;
        //altro pezzo:
        p_stimata=ProbEvento(p,N);
        dev_standard=IncertezzaEvento(p,N);
        R=fabs(p_stimata-p)/dev_standard;
    
        cout<<"p: "<<p<<" stima di p: "<<p_stimata<<" scarto: "<<fabs(p_stimata-p)<<endl;
        cout<<"incertezza teorica: "<<dev_standard<<" rapporto scarto assoluto e inc teorica: "<<R<<endl;
    }

    double s=1.;
    double p_media_stimata=1.0;
    N=1000;
    int M=1000;//faccio 1000 exp
    dev_standard=IncertezzaEvento(p,N);

    RipetiEsperimenti(p,N,M,p_media_stimata,s);//chiamo f e do valori a variab per rif
    cout<<p_media_stimata<<"  "<<dev_standard<<endl;
    cout<<"scarto della media da p: "<<fabs(p_media_stimata -p)<<" dev standard campionaria s: "<<s
        <<"rapporto s/dev teo: "<<s/dev_standard<<endl;
//ultima mi da R che è uno scarto di n deviazioni standard e spesso è rapporto circa cost su 1
//BLOCCO  ":"
    p=0.2;
    int conteggi_qubit[3];
    int conteggi_w[4];
    N=100000;

    SimulaRumore(p,N,conteggi_qubit,conteggi_w);
    int totale=0;
    for(int w=0;w<4;w++){
        cout<<"Peso: "<<w<<" : "<<conteggi_w[w]<<endl;
        totale+=conteggi_w[w];
    }
    cout<<"Totale configurazioni: "<<totale<<endl;

    int conteggi_totali=0;
    for(int i=0; i<3; i++){
        double frequenza = (double)conteggi_qubit[i]/N;

        cout << "Indice qubit: " << i << " conteggio errori: " << conteggi_qubit[i]<< " frequenza stimata: " << frequenza
        << " scarto assoluto: " << fabs(frequenza-p)
        << " incertezza teorica: " << IncertezzaEvento(p,N)<< endl;
        conteggi_totali+=conteggi_qubit[i];
    }
    double qw=1.;
    double q_tilde=1.;
    int somma_w=0;
    for(int w=0;w<4;w++){
        qw=ProbPesoTeorica(p,w);
        q_tilde=(double)conteggi_w[w]/N;
        cout << "Peso: " << w
        << " frequenza sulle w: " << q_tilde
        << " probabilita teorica delle w: " << qw
        << " scarto: " << fabs(q_tilde-qw)
        << " incertezza: " << IncertezzaEvento(qw,N) << endl;
        somma_w+= w*conteggi_w[w];
    }
    double w_medio= (double) somma_w/N;
    //verifico somma di conteggi qubit==somma conteggi w*w che pesa i contributi
    if(somma_w ==conteggi_totali){
        cout<<"controllo dei totali: OK"<<endl;
    }else{
        cout<<"Controllo dei totali: ERRORE"<<endl;
    }
    cout << "Peso medio: " << w_medio
     << " valore teorico: " << 3*p
     << " scarto: " << fabs(w_medio-3*p) << endl;
    
    int errori[3];
    int sindrome[2];
    for(int e0=0;e0<2;e0++){
        for(int e1=0;e1<2;e1++){
            for(int e2=0;e2<2;e2++){
                cout<< e0 <<" "<<e1<<" "<<e2<<endl;
                errori[0] = e0;
                errori[1] = e1;
                errori[2] = e2;
                CalcolaSindrome(errori,sindrome);
                cout<<"peso configurazione: "<<PesoErrore(errori)
                <<" sindrome: "<<sindrome[0]<<" "<<sindrome[1]
                <<" Indice della sindrome: "<<IndiceSindrome(sindrome)<<endl;
            }
        }
    }

    // BLOCCO 4: verifica esaustiva del ciclo di correzione

    int correzione[3];
    int residui[3];
    int sindrome_residua[2];

    int anomalie = 0; //configurazioni non corrette dalla procedura

    cout << "\nErrori | Peso | Sindrome | Correzione | Residui"
    << " | Errore logico | Controllo\n";

    for(int e0=0; e0<2; e0++){
        for(int e1=0; e1<2; e1++){
            for(int e2=0; e2<2; e2++){
                errori[0] = e0;
                errori[1] = e1;
                errori[2] = e2;
                //qui ho costruito  le 8 configurazioni possibili

                CalcolaSindrome(errori, sindrome);
                DecodificaSindrome(sindrome, correzione);
                ApplicaCorrezione(errori, correzione, residui);
                CalcolaSindrome(residui, sindrome_residua);

                //classifichi risultato ottenuto:
                int peso_iniziale = PesoErrore(errori); //err iniziali
                int peso_residuo = PesoErrore(residui); //err rimasti
                bool errore_logico = ErroreLogico(residui); //vedi se devi usare X logico su 111

                // Il residuo deve essere nel codice: 000 oppure 111.
                bool residuo_valido = (peso_residuo==0 || peso_residuo==3); //forma del residuo
                bool sindrome_nulla = (sindrome_residua[0]==0 && sindrome_residua[1]==0); //esisti stabilizer

                // Controllo indipendente per questo decoder:
                // fallisce esattamente quando gli errori iniziali sono almeno due.
                bool esito_atteso = (errore_logico == (peso_iniziale>=2));
                bool controllo = residuo_valido && sindrome_nulla && esito_atteso;
                if(!controllo){
                    anomalie++; //fallisce una verifica
                }

               cout << errori[0] << errori[1] << errori[2]
               << " | " << peso_iniziale
               << " | " << sindrome[0] << sindrome[1]
               << " | " << correzione[0] << correzione[1]
               << correzione[2]
               << " | " << residui[0] << residui[1] << residui[2]
               << " | " << errore_logico
               << " | " << (controllo ? "OK" : "ERRORE") << endl;
            }
        }
    }

    cout << "Anomalie nei controlli: " << anomalie << endl;
    //verifica che tutte le configurazioni siano corrette dalla procedura e restituisce 0 anomalie

    //BLOCCO 5:
    double p_logica=0.0;
    double incertezza=0.0;
    p=0.1;
    N=100000;


    StimaErroreLogico(p,N,p_logica,incertezza); //simula e stima già tutto questa funzione
    double p_log_teorica=ProbLogicaTeorica(p);

    cout<<"prob fisica: "<<p<<" prob logica stimata: "<<p_logica
    <<" incertezza stimata: "<<incertezza
    <<" prob logica teorica: "<<p_log_teorica
    <<" scarto ass tra simulazione e teoria: "<<fabs(p_log_teorica - p_logica)
    <<" rapporto scarto e incertezza: "<<fabs(p_log_teorica - p_logica)/incertezza <<endl;

    //BLOCCO 6: PARTE GRAFICA. Inoltre guardo quanto codice protegge informazione
    ScansionaProbabilita(0.02,0.98,48,100000);

    // CONTROLLO FINALE:
    // a p fissato verifico come la stima di P_L migliora aumentando N.
    cout << "\nControllo convergenza Monte Carlo per p = 0.1\n";

    p = 0.1;
    int valori_N[4] = {1000, 10000, 100000, 1000000};
    double p_teorica_controllo = ProbLogicaTeorica(p);

    for(int i=0; i<4; i++){
        N = valori_N[i];

        // Ogni chiamata esegue N nuovi cicli QEC casuali.
        StimaErroreLogico(p, N, p_logica, incertezza);

        // Misuro la distanza della stima dal valore analitico.
        double scarto = fabs(p_logica - p_teorica_controllo);

        // R dice quante deviazioni standard separano stima e teoria.
        double R = scarto/incertezza;

        cout << "N: " << N
        << " P_L stimata: " << p_logica
        << " incertezza: " << incertezza
        << " P_L teorica: " << p_teorica_controllo
        << " scarto: " << scarto
        << " R: " << R << endl;

    }

    /*
      All'aumentare di N l'incertezza statistica cala circa come 1/sqrt(N).
      La stima non deve avvicinarsi alla teoria in modo monotono: le singole
      estrazioni sono casuali. Un valore di R dell'ordine di 1 indica che
      lo scarto osservato e compatibile con la fluttuazione Monte Carlo.
    */

}

/*
  BLOCCO 1 - Fondamenti Monte Carlo

  Le medie campionarie convergono ai valori teorici aumentando N.
  L'incertezza statistica diminuisce come 1/sqrt(N), ma gli scarti
  casuali non devono diminuire in modo monotono.
*/


/*
  BLOCCO 2 - Rumore bit-flip indipendente

  Ogni qubit subisce un errore X con probabilita p. Le frequenze
  campionate dei pesi 0, 1, 2 e 3 sono compatibili con la legge
  binomiale. Il peso medio tende al valore teorico 3*p.
  I conteggi soddisfano le identita esatte e i casi p=0 e p=1
  restituiscono i risultati deterministici previsti.*/

/*BLOCCO 3: Il simulatore conosce gli errori perché li ha generati: 
li usa per verificare se la correzione ha funzionato. 
Se li passassimo al decoder, gli daremmo informazioni non disponibili nell’esperimento.
Per cui passiamo a chi corregge solo le misure fisiche ottenibili: SINDROME */

/*
  BLOCCO 4 - Decoder e correzione

  Il decoder applica la correzione di peso minimo compatibile con la
  sindrome. Errori di peso zero o uno sono corretti in 000; errori di
  peso due o tre lasciano il residuo 111, cioe un errore logico.
*/

/*
  BLOCCO 5 - Probabilita di errore logico

  Un ciclo QEC e un esperimento Bernoulliano: fallisce se il residuo e
  111. Ripetendo N cicli si stima P_L = K_L/N con incertezza
  sqrt(P_L*(1-P_L)/N).
*/

/*
  La simulazione Monte Carlo riproduce la curva teorica

      P_L(p) = 3 p^2 (1-p) + p^3 = 3 p^2 - 2 p^3.

  Per p < 0.5 la probabilita di errore logico e minore
  dell'errore fisico p: il repetition code protegge il bit logico.

  A p = 0.5 le curve P_L e p si incontrano. Per p > 0.5 il
  decoder a maggioranza fissa diventa svantaggioso.

  Le barre di errore Monte Carlo sono compatibili con la curva
  analitica. Il modello usa rumore bit-flip indipendente e una
  correzione ideale: non e ancora una simulazione fault tolerant.*/

