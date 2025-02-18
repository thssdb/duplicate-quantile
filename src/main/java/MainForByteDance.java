import java.io.IOException;
import java.util.PriorityQueue;

public class MainForByteDance {
    PriorityQueue<Integer> Left,Right;
    void insert(int x){
        if(x<Left.peek()){
            int pop=Left.peek();
            Left.remove(pop);
            Right.add(pop);
            Left.add(x);
        }else{
            int pop=Right.peek();
            Right.remove(pop);
            Left.add(pop);
            Right.add(x);
        }
        if(Left.size()>Right.size()+1){
            int pop=Left.peek();
            Right.add(pop);
            Left.remove(pop);
        }
        if(Right.size()>Left.size()+1){
            int pop=Right.peek();
            Left.add(pop);
            Right.remove(pop);
        }
    }
    int query(){
        if(Left.size()>Right.size())return Left.peek();
        if(Right.size()>Left.size())return Right.peek();
        return (Left.peek()+Right.peek())/2;
    }


    public void main(String[] args) throws IOException {
        Left=new PriorityQueue<>((x,y)->(-Integer.compare(x,y)));
        Right=new PriorityQueue<>((x,y)->(Integer.compare(x,y)));

    }
}
