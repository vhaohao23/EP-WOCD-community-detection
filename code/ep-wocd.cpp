#include<bits/stdc++.h>
using namespace std;

struct Whale{
    vector<int> l;
    vector<int> dk;
    vector<int> lk;
};

int pop=100;
const int T=100;
int N;
int NE;
const double p=1.0;
const double lenP=5.0;

vector<vector<bool>> A;
vector<vector<int>> e;
vector<int> xBest;
vector<int> d;
vector<Whale> whales(pop+1);

random_device rd;   
mt19937 gen(rd());

double modularity(vector<int> dk,vector<int> lk){
    double Q=0;

    for (int i=1;i<=N;i++){
        Q+=double(lk[i])/double(NE)-pow(double(dk[i])/double(2*NE),2.0);
    }

    return Q;
}

void caldklk(vector<int> &l,vector<int> &dk,vector<int> &lk){
    int s=*max_element(l.begin(),l.end());
    dk.assign(N+1,0); lk.assign(N+1,0);
                          

    for (int u=1;u<=N;u++){
        dk[l[u]]+=d[u];
         
        for (int v:e[u])
            if (l[u]==l[v]&&u<v){
                ++lk[l[u]];
            }
    }
}  

void standardization(vector<int> &l,vector<int> &dk,vector<int> &lk){
    map<int,int> mp;
    int cnt=0;
    for (int i=1;i<=N;i++){
        if (mp.find(l[i])==mp.end()){
            ++cnt;
            mp[l[i]]=cnt;
        }
        l[i]=mp[l[i]];
    }

    caldklk(l,dk,lk);
}


void LAR_rand(vector<vector<int>> &a){
    for (int u=1;u<=N;u++){
        if (!e[u].size()) continue;
        uniform_int_distribution<int> disv(0,e[u].size()-1);
        int v=e[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

vector<int> decoding(vector<vector<int>> a){
    bool dd[N+1]={};
    vector<int> l(N+1);
    int cnt=0;

    for (int i=1;i<=N;i++)
        if (!dd[i]){
            ++cnt;
            queue<int> q;
            q.push(i);
            while (!q.empty()){
                int u=q.front();
                q.pop();
                l[u]=cnt;
                for (int v:a[u])
                    if (!dd[v]){
                        dd[v]=true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization(){
    int s;
    for (int p=1;p<=pop;p++){
        vector<vector<int>> a(N+1);
        LAR_rand(a);
        whales[p].l=decoding(a);         
        
        s=*max_element(whales[p].l.begin(),whales[p].l.end());
        whales[p].dk.resize(N+1,0);
        whales[p].lk.resize(N+1,0);                            
        
        for (int u=1;u<=N;u++){
            whales[p].dk[whales[p].l[u]]+=d[u];
             
            for (int v:e[u])
                if (whales[p].l[u]==whales[p].l[v]&&u<v){
                    ++whales[p].lk[whales[p].l[u]];
                }
        }
    }    
}

void movingToPrey(vector<int> &l,vector<int> &dk,vector<int> &lk,double k){
    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        l[pos[i]]=xBest[pos[i]];
    }
    caldklk(l,dk,lk);
}

void randomWalk(vector<int> &l,vector<int> &dk,vector<int> &lk,double k){
    vector<int> ranPop;
    for (int i=1;i<=pop;i++)
        ranPop.push_back(i);
    shuffle(ranPop.begin(),ranPop.end(),gen);   

    vector<vector<int>> P;
    for (int i=0;i<=lenP-1;i++)
        P.push_back(whales[ranPop[i]].l);
    

    vector<int> randNode;
    for (int i=1;i<=N;i++)
        randNode.push_back(i);
    shuffle(randNode.begin(),randNode.end(),gen);


    for (int i=lenP-1;i>=0;i--){
        int kCommma=int(k)/(i+1)+int(int(k)%(i+1)>0);

        while (kCommma--){
            --k;
            l[randNode[k]]=P[i][randNode[k]];
        }
    }
    caldklk(l,dk,lk);
    
}

void encirlingThePrey(vector<int>&l,vector<int> &dk,vector<int> &lk,double r){
    int k=int(r*double(N));

    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++){
        l[pos[i]]=xBest[pos[i]];
    }
    caldklk(l,dk,lk);
}

void mutation(vector<int> &l,vector<int> &dk,vector<int> &lk,double& totalMutation){
    uniform_int_distribution<int> disn(1,N);
    int u=disn(gen);
    
    uniform_real_distribution<double> disType(0.0,1.0);
    double type=disType(gen);

    int S= *max_element(l.begin()+1,l.end());
    ++S;

    if (type<0.5){
        double nodeMuProb=0.5; // probability of mutating each node in the same community
        uniform_real_distribution<double> disProb(0.0,1.0);
        for (int i=1;i<=N;i++)
            if (l[i]==l[u]){
                double randProb=disProb(gen);
                if (randProb<nodeMuProb){
                    l[i]=S;
                }
            }
    }
    else {
        l[u]=S;
        for (int v:e[u]){
            if (l[v]==l[u]){
                l[v]=S;
            }
        }
    }
}

bool isBoundaryNode(vector<int> l,int u){
    for (int v:e[u])
        if (l[u]!=l[v])
            return true;
    return false;
}

void boudaryNodeAdjustment(vector<int> &l,vector<int> &dk,vector<int> &lk){
    vector<int> tmpl;
    vector<int> dktmp;
    vector<int> lktmp;
    int s=*max_element(l.begin(),l.end());
    bool dd[s+1]={};
    
    double inv_4NE2 = 1.0 / (4.0 * double(NE) * double(NE));
    double inv_NE = 1.0 / double(NE);

    for (int i=1;i<=N;i++){
        if (isBoundaryNode(l,i)){
            dd[l[i]]=true;
            for (int neighbor:e[i])
                if (l[i]!=l[neighbor] && !dd[l[neighbor]]){
                    dd[l[neighbor]]=true;
                    
                    tmpl=l;
                    dktmp=dk;
                    lktmp=lk;


                    dk[l[i]]-=d[i];
                    dk[l[neighbor]]+=d[i];
                    for (int v:e[i]){
                        if (l[i]==l[v])
                            --lk[l[i]];
                        if (l[neighbor]==l[v])
                            ++lk[l[neighbor]];
                    }
                    l[i]=l[neighbor];

                    if (modularity(dk,lk)<modularity(dktmp,lktmp)){
                        dk=dktmp;
                        lk=lktmp;
                        l=tmpl;
                    }
                }

            for (int i=1;i<=s;i++)
                dd[i]=false;
        }
    }
}
void EPD(){
    if (whales.size()<10) return;

    vector<pair<double, int>> modularityValues;
    for (int i = 1; i <= pop; i++) {
        double modValue = modularity(whales[i].dk,whales[i].lk);  
        modularityValues.push_back({modValue, i});
    }

    sort(modularityValues.begin(), modularityValues.end());

    vector<Whale> sortedWhales(pop + 1);
    for (int i = 0; i < pop; i++) {
        sortedWhales[i + 1] = whales[modularityValues[i].second];
    }

    whales=sortedWhales;

    double N_nor=pop-(pop/2+1)+1;
    uniform_real_distribution<double> dis(0,1);
    for (int i=pop/2+1;i<=pop;i++){
        double C=1.0-exp(-double(i)/N_nor);
        double rand=dis(gen);
        if (rand<=C){
            whales.erase(whales.begin() + i);
            --pop;
        }
    }
    
}

void updateLocation(vector<int> &l,int t,vector<int> &dk,vector<int> &lk){
    uniform_real_distribution<double> dis(0,1);
    double alpha=dis(gen),beta=dis(gen);
    if (alpha<0.5){
        double k=p*double(t)*double(N)/double(T);

        if (beta<0.5)
            movingToPrey(l,dk,lk,k);
        else
            randomWalk(l,dk,lk,k);
    }
    else {
        uniform_real_distribution<double> disl(-1,1);
        double ll=disl(gen);
        double r=abs(cos(2*M_PI*ll));
        encirlingThePrey(l,dk,lk,r);
    }
}



void EP_WOCD(){
    uniform_real_distribution<double> dis(0,1);
    initialization();
    double ans=0;
    for (int i=1;i<=pop;i++)
        if (modularity(whales[i].dk,whales[i].lk)>ans){
            ans=modularity(whales[i].dk,whales[i].lk);
            xBest=whales[i].l;
        }
    
    int preriodn=5;
    vector<double> mutationProb(pop, dis(gen));
    vector<double> totalMutation(pop, 0);
    vector<double> successMutation(pop, 0);
    
    vector<Whale> whalesTmp=whales;

    for (int t=1;t<=T;t++){
        for (int p=1;p<=pop;p++){
            updateLocation(whales[p].l,t,whales[p].dk,whales[p].lk);
            if (dis(gen)<mutationProb[p]){
                mutation(whales[p].l,whales[p].dk,whales[p].lk,totalMutation[p]);
            }
            standardization(whales[p].l,whales[p].dk,whales[p].lk);
            boudaryNodeAdjustment(whales[p].l,whales[p].dk,whales[p].lk);
            standardization(whales[p].l,whales[p].dk,whales[p].lk);
            
            if (modularity(whales[p].dk,whales[p].lk)>modularity(whalesTmp[p].dk,whalesTmp[p].lk)){
                successMutation[p]++;
            }
        }
        for (int i=1;i<=pop;i++){
            if (modularity(whales[i].dk,whales[i].lk)>ans){
                ans=modularity(whales[i].dk,whales[i].lk);
                xBest=whales[i].l;
            }

            
        }
        EPD();     
        if (t%preriodn==0){
            //do sth
            for (int p=1;p<=pop;p++){
                if (totalMutation[p]){
                    double c=0.9;
                    cout<<successMutation[p]<<" "<<totalMutation[p]<<"\n";
                    if (successMutation[p]/totalMutation[p]>0.2){
                        mutationProb[p]=mutationProb[p]/c;
                    }
                    else if (successMutation[p]/totalMutation[p]<0.2){
                        mutationProb[p]=mutationProb[p]*c;
                    }
                }

                totalMutation[p]=0;
                successMutation[p]=0;
                whalesTmp=whales;
            }
        }
    }    

    cout<<ans<<"\n";
    // for (int i=1;i<=N;i++)
        // cout<<xBest[i]<<" ";
}
int main(){
    clock_t tStart = clock();
    freopen("/home/vhaohao/hao/nckh/dataset-community/karate.txt","r",stdin);


    cin>>N;
    cin>>NE;

    d.resize(N+1);
    A.resize(N + 1, vector<bool>(N + 1, 0));
    e.resize(N+1);
    for (int i=1;i<=NE;i++){
        int u,v;
        cin>>u>>v;
        // u++,v++;
        e[u].push_back(v);
        e[v].push_back(u);
        d[u]++,d[v]++;
        A[u][v]=A[v][u]=true;
    }

    EP_WOCD();

    // printf("\nTime taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}