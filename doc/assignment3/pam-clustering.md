## D: Number of clusters

Using the silhouette score, we can determine the optimal numbers of clusters:

![PAM, optimal number of clusters](../../plots/pam_num_clust.png)

Visually the cluster plots look sparse, no matter the number of clusters or genes. This might be due to 
the pre-normalized Dataset. While using the silhouette method the optimal number of clusters was two, visually four distinct clusters can be determined, which seems to be a better fit.

![PAM k=2 n=10000 clusters](../../plots/pam_k2-n10000.png)

![PAM k=3 n=10000 clusters](../../plots/pam_k3-n10000.png)

![PAM k=4 n=10000 clusters](../../plots/pam_k4-n10000.png)

## Number of genes and how it affects clustering
![PAM k=4 n=10 clusters](../../plots/pam_k4-n10.png)
![PAM k=4 n=100 clusters](../../plots/pam_k4-n100.png)
![PAM k=4 n=1000 clusters](../../plots/pam_k4-n1000.png)
![PAM k=4 n=5000 clusters](../../plots/pam_k4-n5000.png)
![PAM Alluvial plots for n 10,100,1000,5000](../../plots/pam_alluvial.png)


