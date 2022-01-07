#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "hlist.h"

void mkListHash(LINK listHash[],LINK head)
{
	int n,i;
	assert(head != NULL);
	n=listlength(head);
	for(i=0;i<n;++i)
	{
		listHash[i]=head; //
		if(head->next != NULL)
			head=head->next;
	}
}

LINK string_to_list(char s[])
{
	LINK head;
	if(s[0] == '\0')
		return NULL;
	else {
		head = malloc(sizeof(ELEMENT));
		head -> data.c = s[0];
		head -> next = string_to_list(s+1);
		return head;
	}
}

LINK str2list(char s[])
{
	LINK head = NULL,tail,previous;
	int i;
	if(s[0] != '\0') {
		head = malloc(sizeof(ELEMENT));
		head -> data.c = s[0];
		tail = head; previous = head;
		for(i=1; s[i] != '\0'; ++i) {
			tail -> next = malloc(sizeof(ELEMENT));
			tail = tail -> next; tail -> prev = previous;	
			tail -> data.c =s[i]; previous = tail;
		}
		tail -> next = NULL;
	}
	return head;
}

LINK initAtom()
{
	LINK head = NULL,tail,previous;
	head = malloc(sizeof(ELEMENT));
	tail = head; previous = head;
	head -> prev = NULL;
	tail -> next = NULL;
	return head;
}

void addAtom(LINK at)
{
	LINK addingAtom = NULL,tail,tail2;
	addingAtom = malloc(sizeof(ELEMENT));
	tail=at; tail2=at;
	while(tail-> next != NULL){
		tail = tail->next; tail ->prev = tail2;
		tail2= tail;
	}
	tail -> next = addingAtom;
	addingAtom -> prev = tail2;
	addingAtom -> next = NULL;
}

LINK getAtom(LINK at,int nr)
{
	int i; 
	LINK outAt=NULL,tail,tail2;

	tail=at; tail2=at;
	for(i=1;i<nr;i++){
		tail = tail -> next; tail -> prev = tail2;
		tail2= tail;
	}
	return(tail);	

}

LINK getFirst(LINK at)
{
	assert(at!=NULL);
	for ( ; at->prev != NULL ; at = at -> prev)
	return at;
}

void copyAtom(LINK at1,LINK at2)
{
	at1 -> data = at2 ->data;
}


void swapAtoms(LINK at1,LINK at2)
{
	LINK temp=NULL;
	temp = malloc(sizeof(ELEMENT));
	temp -> data = at1 -> data;
	at1  -> data = at2 -> data;
	at2  -> data = temp -> data;
	free(temp);
}

void setFullAtom(LINK at,int nr,char* field, char* atn, char* resn,float x,float y,float z, float mass, float charge,int cgid)
{
	at -> data.cgnr = cgid; 
	at -> data.resnr = nr; 
	strncpy(at -> data.ftyp,field,DIME*2-1);
	strncpy(at -> data.name,atn,DIME); 
	strncpy(at -> data.resi,resn,DIME);
	at -> data.pos[0] = x; at -> data.pos[1] = y; at -> data.pos[2] = z;
	at -> data.name[DIME]='\0'; at -> data.resi[DIME]='\0';
	at -> data.mass_ch[0]=mass; at -> data.mass_ch[1]=charge;
}

void addFullAtom(LINK at,int nr,char* field, char* atn, char* resn,float x,float y,float z, float mass, float charge,int cgid)
{
	LINK addingAtom = NULL,tail,tail2;
	addingAtom = malloc(sizeof(ELEMENT));
	tail=at; tail2=at;
	while(tail-> next != NULL){
		tail = tail->next; tail ->prev = tail2;
		tail2= tail;
	}
	tail -> next = addingAtom;
	addingAtom -> prev = tail2;
	setFullAtom(addingAtom,nr,field,atn,resn,x,y,z,mass,charge,cgid);
	addingAtom -> next = NULL;
}


int iterCount(LINK head)
{
	int cnt = 0;

	for ( ; head != NULL ; head = head -> next)
		++cnt;
	return cnt;
}

void print_listc(LINK head)
{
	if(head == NULL)
		printf("EMPTY");
	else {
		printf("%c --> ", head -> data.c);
		print_listc(head -> next);
	}
}

LINK initFullAtom(int nr,char* field, char* atn, char* resn,float x,float y,float z, float mass, float charge,int cgid)
{
	LINK head = NULL,tail,previous;
	head = malloc(sizeof(ELEMENT));
	tail = head; previous = head;
	setFullAtom(head,nr,field,atn,resn,x,y,z,mass,charge,cgid);
	tail -> next = NULL;

	return head;
}

void printAtom(LINK at)
{
	printf("\n");
	printf("Atom name: %s  Residue: %s  ResID: %d",at -> data.name, at -> data.resi, at -> data.resnr);
	printf("\nCoordinate: %-f %-f %-f",at -> data.pos[0], at -> data.pos[1],at -> data.pos[2]);
	printf("\nForcefield: %s",at -> data.ftyp);
	printf("\nCharge: %f  Mass: %f  Cgnr: %d",at -> data.mass_ch[0],at -> data.mass_ch[1],at -> data.cgnr);
	printf("\nQMMM: %d",at->data.qmmm);
	printf("\n");
}

void printAtomnr(LINK at, int nr)
{
	int i;

	for(i=0;i<nr-1;i++)
		at = at -> next;
	printAtom(at);
}

void printAtomList(LINK head)
{
	if(head == NULL)
		printf("\nEMPTY\n");
	else {
		printAtom(head);
		printf("\n |");
		printf("\n V");
		printAtomList(head -> next);
	}
}

void printlink(LINK p)
{
	assert(p!=NULL);
	printf("%c",p->data.c);
}

void printlinknr(LINK l,int nr)
{
	int i;

	for(i=1;i<nr-1;i++)
		l = l -> next;
	printf("%c",l->data.c);
}

void concatenate(LINK a, LINK b)
{
	assert(a != NULL);
	if(a->next == NULL) {
		a -> next = b;
		b -> prev = a;
	}
	else
		concatenate(a -> next, b);
}

void inAtomnr(LINK l,LINK p,int nr)
{
	// ATOM GETS THE INDEX nr 
	// GETS PLACED BEFORE OLD INDEX

	LINK LOOSE=NULL;
	int i,n,nl;
	n=nr-1; nl=listlength(l);
	
	assert(l!=NULL);
	assert(p!=NULL);

	if(nr <= nl){
		if( nr > 1 ){
			for(i=1;i<n;i++)
				l = l -> next;
			p -> next = l -> next;
			l -> next = p;
			p -> next -> prev = p;
		}
		if(nr == 1){
			p -> next = l -> next;
			l -> next = p;
			p -> next -> prev = p;
			swapAtoms(l,p);
		}

	}else	if(nr == nl+1){
		addAtom(l);
		LOOSE=getAtom(l,nl+1);
		LOOSE -> data  = p -> data;
	}
	else
		printf("\n*********\nCANNOT INSERT OUT OF LIST\n***********\n%d %d",nr,nl);
}

void insert(LINK p1, LINK p2, LINK q)
{
	assert(p1 -> next == p2);
	p1 -> next = q;
	q -> prev = p1;
	q -> next = p2;
	p2 -> prev = q;
}

void deletenr(LINK l, int nr)
{
	int i,inl,n;
	LINK start,fini;
	inl=listlength(l);
	n=nr;

	if(n>1)
		for(i=1; i< n; ++i)
			l = l -> next;

	if(n==1){
		l -> next -> prev=NULL; 
		l -> prev = NULL;
	}
	else if(n==inl){
		l->prev->next=NULL;
	}
	else{
	l-> prev -> next = l->next;
	l-> next -> prev = l->prev;
	}

	free(l);
}

void delAtom(LINK l,int nr)
{
	int i,inl,n;
	LINK start,other;
	inl=listlength(l);
	n=nr;

	if(n==1){
		other = getAtom(l,2);
		swapAtoms(l,other);
		delAtom(l,2);
	}
	if(n>1){
		for(i=1; i< n; ++i)
			l = l -> next;
		if(n==inl){
			l->prev->next=NULL;
		}
		else{
			l-> prev -> next = l->next;
			l-> next -> prev = l->prev;
		}
		free(l);
	}
	
	
}

void delete_list(LINK head)
{
	if(head != NULL) {
		delete_list(head -> next);
		free(head);
	}
}

int listlength(LINK head)
{
	if(head == NULL)
		return 0;
	else
		return (1 + listlength(head -> next));	
}
