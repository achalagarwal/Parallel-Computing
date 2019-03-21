package edu.coursera.concurrent;

import edu.rice.pcdp.Actor;
import edu.rice.pcdp.PCDP;

/**
 * An actor-based implementation of the Sieve of Eratosthenes.
 *
 * TODO Fill in the empty SieveActorActor actor class below and use it from
 * countPrimes to determine the number of primes <= limit.
 */
public final class SieveActor extends Sieve {
    /**
     * {@inheritDoc}
     *
     * TODO Use the SieveActorActor class to calculate the number of primes <=
     * limit in parallel. You might consider how you can model the Sieve of
     * Eratosthenes as a pipeline of actors, each corresponding to a single
     * prime number.
     */

    @Override
    public int countPrimes(final int limit) {
        int count = 1;
        final SieveActorActor actor = new SieveActorActor(2);
        PCDP.finish(()-> {

            for (int i = 3; i <= limit; i++) {
                actor.send(i);
            }

        });
        PCDP.finish(()->
            actor.process(0));
            SieveActorActor t = actor;

            while (t != null) {
                count += t.count;
                t = t.nextActor;
            }

        return count;
        //throw new UnsupportedOperationException();
    }

    /**
     * An actor class that helps implement the Sieve of Eratosthenes in
     * parallel.
     */
    public static final class SieveActorActor extends Actor {
        /**
         * Process a single message sent to this actor.
         *
         * TODO complete this method.
         *
         * @param msg Received message
         */
        SieveActorActor nextActor;
        int value;
        int count = 0;


        SieveActorActor(int v){
            this.value = v;
        }

        @Override
        public void process(final Object msg) {

            Integer i = (Integer)msg;

            if(i == 0){
                if(nextActor!=null) {
                    PCDP.async(()->
                        nextActor.send(i));
                }

                //else
                  //  System.exit(0);
            }


            if(i == value){
                count +=1;
                //System.out.println(value);
                return;
            }

            if(i % value == 0){
                return;
            }
            if(nextActor!=null)
                nextActor.send(msg);
            else{
                nextActor = new SieveActorActor(i);
                PCDP.async(()->
                nextActor.send(i));
            }

        }
    }
}
