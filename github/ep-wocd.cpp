#include<bits/stdc++.h>
using namespace std;
int pop=100;
const int T=100;
int N;
int NE;
const double p=1.0;
const double lenP=5.0;

vector<vector<bool>> A;
vector<vector<int>> x(pop+1);
vector<vector<int>> e;
vector<int> xBest;
vector<int> d;

random_device rd;   
mt19937 gen(rd());

double modularity(vector<int> l){
    int S=0;
    double Q=0;
    for (int label:l)
        S=max(S,label);
    
    vector<int> cs[S+1];

    for (int i=1;i<=N;i++)
        cs[l[i]].push_back(i);

    for (auto c:cs){
        double sumd=0;
        for (int u:c){
            sumd+=d[u];
            for (int v:c)
                if (u<v)
                    Q+=double(A[u][v])/double(NE);
        }
        Q-=pow(sumd/(2*NE),2);
    }
    return Q;
}

void LAR_rand(vector<vector<int>> &a){
    for (int u=1;u<=N;u++){
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
    for (int p=1;p<=pop;p++){
        vector<vector<int>> a(N+1);
        LAR_rand(a);
        x[p]=decoding(a);                                       
    }    
}

void movingToPrey(vector<int> &l,double k){
    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++)
        l[pos[i]]=xBest[pos[i]];
}

void randomWalk(vector<int> &l,double k){
    vector<int> ranPop;
    for (int i=1;i<=pop;i++)
        ranPop.push_back(i);
    shuffle(ranPop.begin(),ranPop.end(),gen);   

    vector<vector<int>> P;//need optimize
    for (int i=0;i<=lenP-1;i++)
        P.push_back(x[ranPop[i]]);
    

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
}

void encirlingThePrey(vector<int>&l,double r){
    int k=int(r*double(N));

    vector<int> pos;
    for (int i=1;i<=N;i++)
            pos.push_back(i);

    shuffle(pos.begin(),pos.end(),gen);

    for (int i=0;i<=int(k)-1;i++)
        l[pos[i]]=xBest[pos[i]];
}

void mutation(vector<int> &l,double u){
    uniform_real_distribution<double> dis(0,1);
    
    int S = *max_element(l.begin(), l.end());
    vector<int> ltmp;
    for (int i=1;i<=N;i++){
        double x=dis(gen);
        if (x<u){
            ltmp=l;
            ++S;
            double y=dis(gen);
            if (y<0.5){                
                l[i]=S;
            }
            else{
                l[i]=S;
                for (int neigbor:e[i])
                    l[neigbor]=S;
            }
            
            if (modularity(l)<modularity(ltmp)){
                --S;
                l=ltmp;
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

void boudaryNodeAdjustment(vector<int> &l){
    vector<int> tmpl;
    for (int i=1;i<=N;i++){
        if (isBoundaryNode(l,i)){
            for (int neighbor:e[i])
                if (l[i]!=l[neighbor]){
                    tmpl=l;

                    l[i]=l[neighbor];

                    if (modularity(l)<modularity(tmpl))
                        l=tmpl;
                }
        }
    }
}
void EPD(){
    if (x.size()<10) return;

    vector<pair<double, int>> modularityValues;
    for (int i = 1; i <= pop; i++) {
        double modValue = modularity(x[i]);
        modularityValues.push_back({modValue, i});
    }

    sort(modularityValues.begin(), modularityValues.end());

    vector<vector<int>> sortedX(pop + 1);
    for (int i = 0; i < pop; i++) {
        sortedX[i + 1] = x[modularityValues[i].second];
    }

    x=sortedX;

    double N_nor=pop-(pop/2+1)+1;
    uniform_real_distribution<double> dis(0,1);
    for (int i=pop/2+1;i<=pop;i++){
        double C=1.0-exp(-double(i)/N_nor);
        double rand=dis(gen);
        if (rand<=C){
            x.erase(x.begin() + i);
            --pop;
        }
    }
    
}

void upadateLocation(vector<int> &l,int t){
    uniform_real_distribution<double> dis(0,1);
    double alpha=dis(gen),beta=dis(gen);
    if (alpha<0.5){
        double k=p*double(t)*double(N)/double(T);

        if (beta<0.5)
            movingToPrey(l,k);
        else
            randomWalk(l,k);
    }
    else {
        uniform_real_distribution<double> disl(-1,1);
        double ll=disl(gen);
        double r=abs(cos(2*M_PI*ll));
        encirlingThePrey(l,r);
    }
}

void EP_WOCD(){
    initialization();
    double ans=0;
    for (int i=1;i<=pop;i++)
        if (modularity(x[i])>ans){
            ans=modularity(x[i]);
            xBest=x[i];
        }
        
    for (int t=1;t<=T;t++){
        for (int p=1;p<=pop;p++){
            upadateLocation(x[p],t);
            mutation(x[p],0.3);
            boudaryNodeAdjustment(x[p]);
            
        }
        for (int i=1;i<=pop;i++)
            if (modularity(x[i])>ans){
                ans=modularity(x[i]);
                xBest=x[i];
            }
        EPD();        
    }    

    cout<<ans<<"\n";
    for (int i=1;i<=N;i++)
        cout<<xBest[i]<<" ";
}
int main(){
    freopen("input.txt","r",stdin);

    cin>>N;
    cin>>NE;

    d.resize(N+1);
    A.resize(N + 1, vector<bool>(N + 1, 0));
    e.resize(N+1);
    for (int i=1;i<=NE;i++){
        int u,v;
        cin>>u>>v;
        e[u].push_back(v);
        e[v].push_back(u);
        d[u]++,d[v]++;
        A[u][v]=A[v][u]=true;
    }
    EP_WOCD();;
}
