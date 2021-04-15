// TODO: modify wasserstein_norm
//
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<algorithm>
#include<set>
#include<list>
#include"metrics.hpp"

using namespace std;


// pos1, pos2 are increasing
double __wasserstein1_norm(double *pos1, double *mass1, int size1, double scale1, double *pos2, double *mass2, int size2, double scale2) {
	double now;

	if((size1==0) && (size2==0)) {
		return 0.;
	} else if(size1==0) {
		now=pos2[0];
	} else if(size2==0) {
		now=pos1[0];
	} else {
		now=min(pos1[0], pos2[0]);
	}
	double metric=0.;
	double F=0.;
	double G=0.;
	int p1=0, p2=0;

	while((p1<size1) && (p2<size2)) {
		if(pos1[p1] <= pos2[p2]) {
			metric += (pos1[p1] - now) * abs(F-G);
			now=pos1[p1];
			F+=mass1[p1]*scale1;
			p1++;
		} else if(pos1[p1] > pos2[p2]) {
			metric += (pos2[p2] - now) * abs(F-G);
			now = pos2[p2];
			G+=mass2[p2]*scale2;
			p2++;
		}
	}
	for( ; p1<size1; p1++) {
		metric += (pos1[p1] - now) * abs(F-G);
		now=pos1[p1];
		F+=mass1[p1]*scale1;
	}
	for( ; p2<size2; p2++) {
		metric += (pos2[p2] - now) * abs(F-G);
		now=pos2[p2];
		G+=mass2[p2]*scale2;
	}

	if( abs(F-G) < eps )
		return metric;
	return HUGE_VAL;
}
double wasserstein1(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	__wasserstein1_norm(pos1,mass1,size1,1.0,pos2,mass2,size2,1.0);
}

double wasserstein1_norm(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	double tmass1=0.;
	double tmass2=0.;
	for(int i=0;i<size1;i++) 
		tmass1+=mass1[i];
	for(int i=0;i<size2;i++) 
		tmass2+=mass2[i];
	
	if(tmass1<tmass2) {
		return abs(tmass1-tmass2) + __wasserstein1_norm(pos1, mass1,size1,tmass2/tmass1,pos2,mass2,size2,1.0);
	} else {
		return abs(tmass1-tmass2) + __wasserstein1_norm(pos2, mass2,size2,tmass1/tmass2,pos1,mass1,size1,1.0);
	}
}

// pos1, pos2 are increasing
double radon(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	double metric=0.;
	int p1=0, p2=0;

	while((p1<size1) && (p2<size2)) {
		if(pos1[p1] < pos2[p2]) {
			metric += mass1[p1];
			p1++;
		} else if(pos1[p1] > pos2[p2]) {
			metric += mass2[p2];
			p2++;
		} else {
			metric += abs(mass1[p1]-mass2[p2]);
			p1++;
			p2++;
		}
	}
	for( ; p1<size1; p1++) {
		metric += mass1[p1];
	}
	for( ; p2<size2; p2++) {
		metric += mass2[p2];
	}
	return metric;
}

double __flat(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	double *fnx = (double*)malloc( (size1+size2+3)*sizeof(double));
	double *fnf = (double*)malloc( (size1+size2+3)*sizeof(double));
	int fnn;

	double prev_pos = -HUGE_VAL;
	double prev_mass=0.;
	double next_pos;
	double next_mass;
	double most_left;

	int p1=0, p2=0;
	
	while( (p1<size1) || (p2<size2) ) {
		if( (p1<size1)&&(p2<size2) ) {
			if(pos1[p1] < pos2[p2]) {
				next_pos=pos1[p1];
				next_mass=mass1[p1];
				p1++;
			} else if (pos1[p1] == pos2[p2] ) {
				next_pos=pos1[p1];
				next_mass=mass1[p1]-mass2[p2];
				p1++;
				p2++;
			} else {
				next_pos=pos2[p2];
				next_mass=-mass2[p2];
				p2++;

			}
		} else if (p1<size1) {
			next_pos=pos1[p1];
			next_mass=mass1[p1];
			p1++;
		} else {
			next_pos=pos2[p2];
			next_mass=-mass2[p2];
			p2++;
		}
		//printf("p=(%d,%d) now=(%lf,%lf) prev=(%lf,%lf)\n", p1,p2,next_pos,next_mass,prev_pos,prev_mass);
		if(prev_pos == -HUGE_VAL) {
			fnx[0] = -1;
			fnx[1] = 1;
			fnf[0] = next_mass;
			fnf[1] = -HUGE_VAL;
			most_left = -next_mass;
			fnn=2;
		} else {
			double d = next_pos - prev_pos;
			int w,k;
			double val = most_left;

			w=0;				// where to write new sequence
			for(k=0; k<fnn; k++) {
			//	printf("k=%d a=%lf d=%lf fnn=%d\n",k,fnf[k],d,fnn);

				if(fnx[k]-d >= -1.0) {
					if(w==0) {	// havent been in [-1,1] yet
						val += (-1 - (fnx[k-1]-d)) * fnf[k-1];
						most_left=val;
						fnx[w] = -1;
						fnf[w] = fnf[k-1];

						w++;
					} 
					if(fnf[k]>=0.) {
						fnx[w]=fnx[k]-d;
						fnf[w]=fnf[k];
						w++;
					}
				} else {
					if(k>0)
						val += (fnx[k]-fnx[k-1]) * fnf[k-1];
				}
				if(fnf[k] < 0)
					break;
				
			}
		/*	
			printf("Got here iwth w=%d, k=%d, fnn=%d val=%lf mostleft=%lf\n", w,k,fnn,val,most_left);
			for(int l=0;l<fnn;l++) {
				printf("\t\tx=%lf a=%lf\n", fnx[l],fnf[l]);
			}
			*/
			

			if(w==0){	// havent been in [-1,1] -> max will propagate to -1
				most_left=val;
			}
			for(int j=fnn; j>k; j--) { // might cause SEGV if two deltas are in the same place
				fnf[j]=fnf[j-1];
				fnx[j]=fnx[j-1]+d;
			}
			fnn++;
			fnf[k]=0.;
			fnx[k]=max(-1.0,fnx[k]-d);
			most_left += next_mass*fnx[0];
			for(int j=k; j<fnn; j++) {
				fnf[w]=fnf[j];
				fnx[w]=fnx[j];
				if(fnx[w] >= 1.0) {
					fnx[w]=1.0;
					fnf[w]=-HUGE_VAL;
					w++;
					break;
				}
				w++;
			}
			fnn=w;
			for(int i=0;i<fnn;i++) {
				fnf[i]+=next_mass;
			}
		}

		prev_pos=next_pos;
		prev_mass=next_mass;
		
		double val=most_left;
		int i=0;
		/*
		printf("<plot>\n");
		printf("%lf %lf %lf\n",fnx[0],most_left,fnf[0]);
		do
		{
			val += (fnx[i+1] - fnx[i]) * fnf[i];
			printf("%lf %lf %lf\n", fnx[i+1], val,fnf[i+1]);
			i++;
		} while ( (fnf[i] != -HUGE_VAL) ) ;
		printf("e\n</plot>\n");
		*/
		
	}
	double val=most_left;
	int i=0;
	while( (fnf[i] >= 0) && (fnx[i]<-1.)) {
		val += (fnx[i+1] - fnx[i]) * fnf[i];
		i++;
	}
	return val;
}

double __flat_nlogn(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	set<fentry> f;

	bool first_roll=true;
	double d;
	double prev_pos = -HUGE_VAL;
	double prev_mass=0.;
	double next_pos;
	double next_mass;

	double left_value;
	double a_mod=0.;

	int p1=0, p2=0;
	
	while( (p1<size1) || (p2<size2) ) {
		if( (p1<size1)&&(p2<size2) ) {
			if(pos1[p1] < pos2[p2]) {
				next_pos=pos1[p1];
				next_mass=mass1[p1];
				p1++;
			} else if (pos1[p1] == pos2[p2] ) {
				next_pos=pos1[p1];
				next_mass=mass1[p1]-mass2[p2];
				p1++;
				p2++;
			} else {
				next_pos=pos2[p2];
				next_mass=-mass2[p2];
				p2++;

			}
		} else if (p1<size1) {
			next_pos=pos1[p1];
			next_mass=mass1[p1];
			p1++;
		} else {
			next_pos=pos2[p2];
			next_mass=-mass2[p2];
			p2++;
		}
		if(first_roll) {
			f.insert(fentry(-1., next_mass));
			f.insert(fentry(2,-HUGE_VAL));
			left_value=-next_mass;
			first_roll=false;
			prev_pos=next_pos;
			prev_mass=next_mass;
			continue;
		}
		d = next_pos - prev_pos;
		/*
		printf("BEFORE(%lf,%lf)\nleft=%lf amod=%lf d=%lf\n", next_pos,next_mass,left_value,a_mod,d);
		for(set<fentry>::reverse_iterator rit=f.rbegin(); rit!=f.rend(); rit++) {
			printf("x=%lf a=%lf\n",rit->x,rit->a);
		}
		*/

		set<fentry>::iterator fdec = --f.lower_bound(fentry(-a_mod)); // pierwszy malejacy
		set<fentry>::iterator fdec_n=fdec;
		set<fentry>::iterator fdec_p=fdec;

		fdec_n++;
		fdec_p--;

		// FROM: (fdec_n) fdec (fdec_p)
		// TO:   (fdec_n)-d (0)-d fdec+d (fdec_p)+d
		// FROM: fdec (fdec_p)
		// TO:   (-1) fdec+d (fdec+1)+d
		if( (fdec_n) == f.end() ) {
			if((fdec->a) != (-a_mod)) {
				*(const_cast<double*>(&(fdec->x))) += 1+2*d; // fdec->x += 1+2*d;
				f.insert(fentry(-1.-d,(-a_mod)));
			} else {
				*(const_cast<double*>(&(fdec_p->x))) += 2*d;
			}
		} else {
			*(const_cast<double*>(&((--f.end())->x))) -=d;
			if( (fdec->a) != (-a_mod) ) {
				f.insert(fentry(fdec->x,(-a_mod)));
				*(const_cast<double*>(&(fdec->x))) = 2*d;
			} else {
				*(const_cast<double*>(&(fdec_p->x))) += 2*d;
			}
		}
		/*
		printf("MIDDLE(%lf,%lf)\nleft=%lf amod=%lf\n", next_pos,next_mass,left_value,a_mod);
		for(set<fentry>::reverse_iterator rit=f.rbegin(); rit!=f.rend(); rit++) {
			printf("x=%lf a=%lf\n",rit->x,rit->a);
		}
		*/

		double pos=1+d;;
		set<fentry>::iterator it;
		for(it=f.begin(); it!=f.end(); it++) {
			pos -= it->x;
			if(pos<1.) {
				break;
			}
		}
		f.erase(f.begin(), ++it);
		f.insert(fentry(1.-pos,-HUGE_VAL));
		
		/*
		printf("FIXEDEND(%lf,%lf)\nleft=%lf amod=%lf\n", next_pos,next_mass,left_value,a_mod);
		for(set<fentry>::reverse_iterator rit=f.rbegin(); rit!=f.rend(); rit++) {
			printf("x=%lf a=%lf\n",rit->x,rit->a);
		}
		*/


		double val=left_value;
		set<fentry>::reverse_iterator rit;
		double ppos=-1-d;
		double pdir=0;
		bool quit=false;
		pos=0.;
		it=(--f.end());
		while(!quit) {
			pos+=it->x;
			if(pos>-1) {
				left_value=val+pdir*(-1-ppos);
				*(const_cast<double*>(&(it->x))) =(pos+1);
				it++;
				*(const_cast<double*>(&(it->x))) = -1;
				break;

			}
			val+=pdir*(pos-ppos);
			ppos=pos;
			pdir=it->a + a_mod;
			if(it==f.begin())
				quit=true;
			else
				it--;
		} 
//		printf("   deleting (%lf %lf).next\n",it->x, it->a);
		f.erase(++it,f.end());
		
		a_mod+=next_mass;
		left_value += (-next_mass);

/*
		printf("ALLDONE: left=%lf amod=%lf\n", left_value,a_mod);
		for(set<fentry>::reverse_iterator rit=f.rbegin(); rit!=f.rend(); rit++) {
			printf("x=%lf a=%lf\n",rit->x,rit->a);
		}
		*/
		
		prev_pos=next_pos;
		prev_mass=next_mass;
	}

	double val=left_value;
	double ppos=-1;
	double pos=0;
	double pdir=0;
	for(set<fentry>::reverse_iterator rit=f.rbegin(); rit!=f.rend(); rit++) {
		pos+=rit->x;
		val+=pdir*(pos-ppos);
		ppos=pos;
		pdir=rit->a + a_mod;
		if(pdir<0)
			return val;
	}
	return HUGE_VAL;
}

double __flat_n(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {

	double tmass=0.;
	for(int i=0;i<size1; i++)
		tmass+=mass1[i];
	for(int i=0;i<size2; i++)
		tmass-=mass2[i];

	int size = size1+size2;
	double *pos = (double*) malloc(size*sizeof(double));
	double *mass = (double*) malloc(size*sizeof(double));
	double *mod = (double*) calloc(size,sizeof(double));
	list<int> neigh;
	double now;

	bool posmod = (tmass>0); // true if 1 larger than 2
	int p1=0, p2=0;

	if((size1==0) && (size2==0)) {
		return 0.;
	} else if(size1==0) {
		now=pos2[0];
	} else if(size2==0) {
		now=pos1[0];
	} else {
		now=min(pos1[0], pos2[0]);
	}

	while((p1<size1) && (p2<size2)) {
		if(pos1[p1] <= pos2[p2]) {
			pos[p1+p2] = pos1[p1];
			mass[p1+p2] = (posmod)?-mass1[p1]:mass1[p1];
			p1++;
		} else if(pos1[p1] > pos2[p2]) {
			pos[p1+p2] = pos2[p2];
			mass[p1+p2] = (posmod)?mass2[p2]:-mass2[p2];
			p2++;
		}
	}
	for( ; p1<size1; p1++) {
		pos[p1+p2] = pos1[p1];
		mass[p1+p2] = (posmod)?-mass1[p1]:mass1[p1];
	}
	for( ; p2<size2; p2++) {
		pos[p1+p2] = pos2[p2];
		mass[p1+p2] = (posmod)?mass2[p2]:-mass2[p2];  
	}

	for(int p=0; p<size; p++) {
		if(mass[p] < 0)
			neigh.push_back(p);
	}

	list<int>::iterator it=neigh.begin();
	list<int>::iterator it2;
	double dp,dm,m;
	for(int p=0; p<size; p++) {

		while((it!=neigh.end()) && (pos[*it] < pos[p]))
			it++;

		while(mass[p]>0.) {
			if(neigh.begin() == neigh.end())
				break;
			it2 = it;
			it2--;

			dp = (it==neigh.end()) ? HUGE_VAL : abs(pos[*it]-pos[p]);
			dm = (it2==neigh.end()) ? HUGE_VAL : abs(pos[*it2]-pos[p]);
			
			if(dp<dm) {
				m = min(mass[p], -mass[*it]);
				printf("moving %lf from %lf to %lf\n", m, pos[p], pos[*it]);
				mod[p] += m;
				mod[*it] -= m; 
				mass[p]-=m;
				mass[*it]+=m;

				if(mass[*it]==0.)
					neigh.erase(it);
				it=it2;
				it++;

			} else {
				m = min(mass[p], -mass[*it2]);
				printf("moving2 %lf from %lf to %lf\n", m, pos[p], pos[*it2]);
				mod[p] += m;
				mod[*it2] -= m;
				mass[p]-=m;
				mass[*it2]+=m;
				if(mass[*it2]==0.)
					neigh.erase(it2);
			}


		}

	}
	/*
	printf("mod=\n");
	for(int p=0;p<size;p++)
		printf("%lf ", mod[p]);
	printf("\n");
	*/

	now=pos[0];
	double diff=0.;
	double metric=0.;

	for(int p=0;p<size; p++) {
		metric+=(pos[p]-now)*abs(diff);
//		printf("%lf@      %lf ---> %lf mass %lf\n", mass[p],now, pos[p],diff);
		diff+=mod[p];
		now=pos[p];
	}

	return metric + abs(tmass);
}

double flat(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	return max( __flat(pos1,mass1,size1,pos2,mass2,size2), __flat(pos2,mass2,size2,pos1,mass1,size1) );
}

double flat_nlogn(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	return max( __flat_nlogn(pos1,mass1,size1,pos2,mass2,size2), __flat_nlogn(pos2,mass2,size2,pos1,mass1,size1) );
}
double flat_n(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
//	return max( __flat_n(pos1,mass1,size1,pos2,mass2,size2), __flat_n(pos2,mass2,size2,pos1,mass1,size1) );
	return __flat_n(pos1,mass1,size1,pos2,mass2,size2);
}

double wasserstein1_restricted(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2) {
	double distance=0.;
	double partialFront=0.;
	double partialBack=0.;
	double now;
	
	if((size1==0) && (size2==0)) {
		return 0.;
	} else if(size1==0) {
		now=pos2[0];
	} else if(size2==0) {
		now=pos1[0];
	} else {
		now=min(pos1[0], pos2[0]);
	}
	
	int p1=0, p2=0;
	int toofar=0;

//	printf("here now=%lf\n",now);

	while((p1<size1) && (p2<size2)) {
		if(pos1[p1] <= pos2[p2]) {
			if(pos1[p1]>0.) {
				distance += (0-now)*abs(partialFront);
				toofar=1;
				break;
			}
			distance += (pos1[p1] - now) * abs(partialFront);
			now=pos1[p1];
			partialFront+=mass1[p1];
			p1++;
		} else /* if(pos1[p1] > pos2[p2]) */ {
			if(pos2[p2]>0.) {
				distance += (0-now)*abs(partialFront);
				toofar=1;
				break;
			}
			distance += (pos2[p2] - now) * abs(partialFront);
			now = pos2[p2];
			partialFront-=mass2[p2];
			p2++;
		}
	}
	if(!toofar) {
		for( ; p1<size1; p1++) {
			if(pos1[p1]>0.) {
				distance += (0-now)*abs(partialFront);
				toofar=1;
				break;
			}
			distance += (pos1[p1] - now) * abs(partialFront);
			now=pos1[p1];
			partialFront+=mass1[p1];
		}
		for( ; p2<size2; p2++) {
			if(pos2[p2]>0.) {
				distance += (0-now)*abs(partialFront);
				toofar=1;
				break;
			}
			distance += (pos2[p2] - now) * abs(partialFront);
			now=pos2[p2];
			partialFront-=mass2[p2];
		}
	}
	if(!toofar) {
		return distance + abs(partialFront);
	}

	// the next position is >= 0 - first sum is computed
	
	if(size1==0) {
		now=pos2[size2-1];
	} else if(size2==0) {
		now=pos1[size1-1];
	} else {
		now=max(pos1[size1-1], pos2[size2-1]);
	}
	
	int r1=size1-1;
	int r2=size2-1;
//	printf("now =%lf; p1=%d p2=%d r1=%d r2=%d partialFront=%lf distance=%lf\n",now,p1,p2,r1,r2,partialFront,distance);

	while( (r1>=p1) && (r2>=p2) ) {
		if(pos1[r1] >= pos2[r2]) {
			distance += (now-pos1[r1])*abs(partialBack);
			now=pos1[r1];
			partialBack+=mass1[r1];
			r1--;
		} else {
			distance += (now-pos2[r2])*abs(partialBack);
			now=pos2[r2];
			partialBack-=mass2[r2];
			r2--;
		}

	}
	for( ;r1>=p1;r1--) {
		distance += (now-pos1[r1])*abs(partialBack);
		now=pos1[r1];
		partialBack+=mass1[r1];
	}
	for( ;r2>=p2;r2--) {
		distance += (now-pos2[r2])*abs(partialBack);
		now=pos2[r2];
		partialBack-=mass2[r2];
	}
//	printf("and now=%lf distance =%lf; partialback=%lf\n",now,distance,partialBack);


	distance += now*abs(partialBack);

	return distance+abs(partialBack+partialFront);

}

