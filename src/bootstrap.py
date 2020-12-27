# function calculate the jaccard coefficient of two vectors, return values
def jaccard(labels1,labels2):
    n11 = n10 = n01 = 0
    n = len(labels1)
    # TODO: Throw exception if len(labels1) != len(labels2)
    for i, j in itertools.combinations(range(n), 2):
        comembership1 = labels1[i] == labels1[j]
        comembership2 = labels2[i] == labels2[j]
        if comembership1 and comembership2:
            n11 += 1
        elif comembership1 and not comembership2:
            n10 += 1
        elif not comembership1 and comembership2:
            n01 += 1
    return float(n11) / (n11 + n10 + n01)

# function that resample a given dataset, return new dataset
def resample(df):
    df_new = df
    allelements = []
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            allelements.append(df.iat[i,j])
    for x in range(len(df.index)):
        for y in range(len(df.columns)):
            df_new.iat[x,y] = random.choice(allelements)
    return (df_new)

# Hi Alex, I think maybe returning a list of original clusters which dissovled too often 
# is more reasonable than returning a data structure according to your description of this function? 
def hierarchical_bootstrap_validation(df,k, method='average',metric='sqeuclidean',n_bootstraps=100):
    
        x = shc.linkage(df,method=method,metric=metric)
        # the list of original clusters 
        xx = fcluster(x, k, criterion='maxclust')
        # record the list of list1
        list2 = []
        for i in range(n_bootstraps):
            #new_df = resample(df)
            new_df = df.sample(len(df),replace=True)
            y = shc.linkage(new_df,method=method,metric=metric)
            # the list of new clusters
            yy = fcluster(y, k, criterion='maxclust')
            # record the max jaccard coefficient of each original clusters 

            list1 = []
            for j in np.unique(xx):
                max_Jc = 0
                for h in np.unique(yy):
                    Jc = jaccard(xx==j,yy==h)
                    if Jc >= max_Jc:
                        max_Jc = Jc
                list1.append(max_Jc)
            list2.append(list1)
            
        # record the original cluster that dissolved too often
        list3 = []
        for u in range(len(np.unique(xx))):
            z = 0
            for v in range(n_bootstraps):
                if list2[v][u] < 0.5:
                    z += 1
            # justify if dissovled too often
            if z>= (len(xx)/2):
                list3.append(xx[u])
        
        return (list2,list3)